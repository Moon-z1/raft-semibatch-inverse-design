function run_extension
% 生成独立目标后，两类代理使用相同训练数据进行搜索。
    % 八目标扩展沿用冻结设置；方案在首次模拟前保存。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    old=fullfile(fileparts(root),'raft_trajectory_pilot','outputs');
    if numel(dir(fullfile(od,'run_reachable_*_gp.mat')))==48 && exist(fullfile(od,'run_elapsed.mat'),'file')
        fprintf('All 48 runs already exist. Original search timing and records preserved. Use verify_revision to audit.\n');
        return;
    end
    maxNumCompThreads(1);
    diary(fullfile(od,'extension.log'));
    timer=tic;
    load(fullfile(old,'protocol_data.mat'),'test','c','grid','w');
    seeds=[101 202 303];
    training=cell(1,3);
    avoid=test.X;
    for si=1:3
        z=load(fullfile(old,sprintf('seed_%d_compact.mat',seeds(si))),'res');
        training{si}=z.res.training;
        assert(size(training{si}.X,1)==160);
        avoid=[avoid;training{si}.X];
    end
    z=load(fullfile(old,'hidden_target_truth.mat'),'truth');
    avoid=[avoid;vertcat(z.truth.x)];
    cf=c;
    cf.RelTol=2e-9;
    cf.AbsTol=1e-13;
    cf.output_dt=5;
    cf.max_step=30;
    if ~exist(fullfile(od,'protocol_data.mat'),'file')
        rng(2026101901,'twister');
        targets=struct([]);
        truth=struct([]);
        logrows=struct([]);
        attempts=0;
        odecalls=0;
        tg=tic;
        while numel(targets)<8
            attempts=attempts+1;
            assert(attempts<=10000);
            x=.6*(2*rand(1,4)-1);
            [~,~,ok]=feed_decode(x,c);
            distance=min(pdist2(x,avoid));
            reason='';
            secs=0;
            accept=false;
            if ~ok
                reason='pump_invalid';
            elseif distance<=1e-6
                reason='control_overlap';
            else
                odecalls=odecalls+1;
                tr=tic;
                try
                    rec=observe(x,cf,grid);
                    secs=toc(tr);
                    if rec.feasible && all(isfinite(rec.F))
                        accept=true;
                        reason='accepted';
                    else
                        reason='process_or_grid_invalid';
                    end
                catch ex
                    secs=toc(tr);
                    reason=['exception:' ex.identifier ':' ex.message];end
                end
                logrows=[logrows,struct('proposal',attempts,'accepted',accept,'reason',reason,'seconds',secs,'min_distance',distance,'x1',x(1),'x2',x(2),'x3',x(3),'x4',x(4))];
                if accept
                    k=numel(targets)+1;
                    targets(k).id=sprintf('reachable_%02d',k+4);
                    targets(k).kind='new_reachable';
                    targets(k).F=rec.F;
                    truth(k).id=targets(k).id;
                    truth(k).x=x;
                    truth(k).record=rec;
                    avoid=[avoid;x];
                end
            end
            target_seconds=toc(tg);
            save(fullfile(od,'hidden_target_truth.mat'),'truth','logrows','attempts','odecalls','target_seconds');
            writetable(struct2table(logrows),fullfile(od,'target_generation_log.csv'));
            save(fullfile(od,'protocol_data.mat'),'targets','c','grid','w');
        else
            load(fullfile(od,'protocol_data.mat'),'targets','c','grid','w');
        end
        % 仅传递目标曲线，不向优化器传递其生成控制。
        clear truth avoid x z rec
        shared_cost=struct([]);
        for si=1:3
            seed=seeds(si);
            D=training{si};
            t=tic;
            shared=fit_models(D,[],true);
            sharedsecs=toc(t);
            shared_cost=[shared_cost,struct('seed',seed,'shared_fit_seconds',sharedsecs,'training_calls_reused',160,'trajectory_rank',size(shared.B,2),'explained_variance',shared.explained)];
            for ti=1:8
                methods={'scalar_gp','trajectory_gp'};
                if mod(ti+si,2)==0
                    methods=fliplr(methods);
                end
                for mi=1:2
                    method=methods{mi};
                    file=fullfile(od,sprintf('run_%s_%d_%s.mat',targets(ti).id,seed,method));
                    if exist(file,'file')
                        fprintf('RESUME preserved %s\n',file);
                        continue;
                    end
                    try
                        run=search_one(D,targets(ti),seed,ti+10,method,shared,c,grid,w);
                        save(file,'run','-v7.3');
                        fprintf('%s seed=%d %s: %.8f -> %.8f, %.3fs\n',targets(ti).id,seed,method,run.baseline.RMSE,run.RMSE,run.total_seconds);
                    catch ex
                        failure=struct('target',targets(ti).id,'seed',seed,'method',method,'identifier',ex.identifier,'message',ex.message,'report',getReport(ex,'extended'));
                        save([file '.failure.mat'],'failure');
                        fprintf(2,'%s\n',failure.report);
                        rethrow(ex);
                    end
                end
            end
        end
        elapsed=toc(timer);
        save(fullfile(od,'run_elapsed.mat'),'elapsed','shared_cost');
        writetable(struct2table(shared_cost),fullfile(od,'shared_fit_costs.csv'));
        fprintf('EXTENSION COMPLETE, %.3f seconds; 48 runs, 1440 search evaluations.\n',elapsed);
        diary off
    end
    function run=search_one(D,target,seed,target_index,name,shared,c,grid,w)
        [base,ib]=best_score(D,target.F,w);
        baseline=struct('x',D.X(ib,:),'RMSE',base,'record',get_record(D,ib));
        t=tic;
        DD=D;
        fit_seconds=0;
        prediction_seconds=0;
        isTrajectory=strcmp(name,'trajectory_gp');
        for round=1:3
            tfit=tic;
            if isTrajectory && round==1
                model=shared;
            else
                model=fit_models(DD,target.F,isTrajectory);
            end
            fit_seconds=fit_seconds+toc(tfit);
            [~,bestidx]=best_score(DD,target.F,w);
            rng(seed+10000*target_index+100*round,'twister');
            pool=sample_controls(2048,c);
            local=[DD.X(bestidx,:)+.15*randn(768,4);DD.X(bestidx,:)+.04*randn(768,4)];
            local=max(-1,min(1,local));
            keep=false(size(local,1),1);
            for jj=1:size(local,1)
                [~,~,keep(jj)]=feed_decode(local(jj,:),c);
            end
            pool=unique([pool;local(keep,:)],'rows','stable');
            pool=pool(min(pdist2(pool,DD.X),[],2)>1e-6,:);
            tp=tic;
            [pred,pm]=predict_models(model,pool,isTrajectory);
            if isTrajectory
                vals=sqrt(((pred-target.F).^2)*w);
            else
                vals=pred;
            end
            violation=sum(max(0,.8-pm),2);
            orderkey=vals;
            orderkey(violation>0)=1+10*violation(violation>0);
            [~,order]=sort(orderkey);
            chosen=[];
            for jj=order(:).'
                if isempty(chosen) || min(pdist2(pool(jj,:),pool(chosen,:)))>.025
                    chosen(end+1)=jj;
                end
                if numel(chosen)==10
                    break;
                end
            end
            assert(numel(chosen)==10);
            prediction_seconds=prediction_seconds+toc(tp);
            new=make_data(pool(chosen,:),c,grid);
            DD=append_data(DD,new);
        end
        total_seconds=toc(t);
        assert(size(DD.X,1)==190);
        [score,ix]=best_score(DD,target.F,w);
        hist=cummin(score_data(DD,target.F,w));
        run=struct('target',target,'method',name,'seed',seed,'target_rng_index',target_index,'baseline',baseline,'best_x',DD.X(ix,:),'best_record',get_record(DD,ix),'RMSE',score,'best10',hist(170),'best20',hist(180),'best30',hist(190),'fit_seconds',fit_seconds,'prediction_seconds',prediction_seconds,'total_seconds',total_seconds,'extra_ode_seconds',sum(DD.seconds(161:end)));
        fields=fieldnames(DD);
        extra=DD;
        for j=1:numel(fields)
            f=fields{j};
            extra.(f)=DD.(f)(161:end,:);
        end
        run.extra=extra;
        run.initial_X=D.X;
        run.initial_F=D.F;
        run.initial_M=D.M;
    end
    function D=make_data(X,c,grid)
        n=size(X,1);
        D=struct('X',X,'F',nan(n,numel(grid)),'M',nan(n,2),'feasible',false(n,1),'seconds',zeros(n,1),'balance',nan(n,1),'pump_ok',false(n,1));
        for j=1:n
            r=observe(X(j,:),c,grid);
            D.F(j,:)=r.F;
            D.M(j,:)=r.M;
            D.feasible(j)=r.feasible;
            D.seconds(j)=r.seconds;
            D.balance(j)=r.balance;
            D.pump_ok(j)=r.pump_ok;
        end
    end
    function r=observe(x,c,grid)
        t=tic;
        [rates,~,ok,~,cc]=feed_decode(x,c);
        r=struct('F',nan(1,numel(grid)),'M',[NaN NaN],'feasible',false,'seconds',0,'balance',NaN,'pump_ok',ok);
        if ok
            o=simulate_reactor(cc,rates);
            r.M=[o.Xmol(end),o.Xmass(end)];
            r.feasible=all(r.M>=[c.min_Xmol,c.min_Xmass]-1e-8) && max(o.V)<=c.max_volume+1e-8;
            r.balance=max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]);
            if o.Xmol(end)>=grid(end)
                valid=o.Xmol>1e-8;
                [u,ix]=unique(o.Xmol(valid));
                f=o.Finstant(valid);
                r.F=interp1(u,f(ix),grid,'linear');
            end
        end
        r.seconds=toc(t);
    end
    function D=append_data(D,A)
        fields=fieldnames(D);
        for j=1:numel(fields)
            f=fields{j};
            D.(f)=[D.(f);A.(f)];
        end
    end
    function r=get_record(D,j)
        r=struct('F',D.F(j,:),'M',D.M(j,:),'feasible',D.feasible(j),'balance',D.balance(j));
    end
    function v=score_data(D,target,w)
        v=sqrt(((D.F-target).^2)*w);
        v(~D.feasible | ~isfinite(v))=Inf;
    end
    function [b,j]=best_score(D,target,w)
        [b,j]=min(score_data(D,target,w));
        assert(isfinite(b));
    end
    function model=fit_models(D,target,trajectory)
        valid=all(isfinite(D.M),2);
        model.mol=gpfit(D.X(valid,:),D.M(valid,1));
        model.mass=gpfit(D.X(valid,:),D.M(valid,2));
        if trajectory
            good=all(isfinite(D.F),2);
            F=min(1-1e-6,max(1e-6,D.F(good,:)));
            % 把(0,1)内的组成映射到实数域；预测后用逆变换还原。
        Z=log(F./(1-F));
            model.center=mean(Z,1);
            Z=Z-model.center;
            % SVD把整条轨迹压缩为少数模式；每个模式系数分别拟合GP。
        [~,S,V]=svd(Z,'econ');
            energy=cumsum(diag(S).^2)/sum(diag(S).^2);
            q=min(12,find(energy>=.9995,1));
            model.B=V(:,1:q);
            model.explained=energy(q);
            % 每行是一条轨迹的模式系数，B的每列是一个轨迹模式。
        coeff=Z*model.B;
            recon=1./(1+exp(-(coeff*model.B.'+model.center)));
            model.reconstruction=sqrt(mean((recon-F).^2,2));
            model.gps=cell(1,q);
            for j=1:q
                model.gps{j}=gpfit(D.X(good,:),coeff(:,j));
            end
        else
            good=D.feasible & all(isfinite(D.F),2);
            w=ones(size(D.F,2),1)/(size(D.F,2)-1);
            w([1,end])=w([1,end])/2;
            y=sqrt(((D.F(good,:)-target).^2)*w);
            model.scalar=gpfit(D.X(good,:),log(max(y,1e-7)));
        end
    end
    function g=gpfit(X,y)
        [X,ia]=unique(X,'rows','stable');
        y=y(ia);
        g.mu=mean(y);
        g.sd=max(std(y),1e-8);
        g.net=fitrgp(X,(y-g.mu)/g.sd,'KernelFunction','ardmatern52','Standardize',true,'Sigma',.002,'ConstantSigma',true,'SigmaLowerBound',1e-6,'FitMethod','exact','PredictMethod','exact');
    end
    function y=gppredict(g,X)
        y=predict(g.net,X)*g.sd+g.mu;
    end
    function [F,M]=predict_models(model,X,trajectory)
        M=[gppredict(model.mol,X),gppredict(model.mass,X)];
        if trajectory
            A=zeros(size(X,1),numel(model.gps));
            for j=1:numel(model.gps)
                A(:,j)=gppredict(model.gps{j},X);
            end
            F=1./(1+exp(-(A*model.B.'+model.center)));
        else
            F=exp(gppredict(model.scalar,X));
        end
    end
