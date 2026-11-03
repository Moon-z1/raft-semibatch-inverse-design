% 生成或读取轨迹数据库，比较标量与轨迹代理；已有记录按原规则复用。
function run_pilot
% 先建立共同轨迹数据库，再比较各方法推荐的进料方案。
    % 本阶段使用独立保存的核心函数，避免与其他阶段同名函数混用。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    timer=tic;
    maxNumCompThreads(1);
    c=study_config;
    seeds = [101,202,303];
    filesExist = true;
    for k = 1:numel(seeds)
        fileName = sprintf('seed_%d.mat',seeds(k));
        if exist(fullfile(od,fileName),'file') == 0
            filesExist = false;
        end
    end
    if filesExist
        fprintf('All seed results already exist. Original timing preserved. Run analyze_pilot to rebuild summaries.\n');
        return;
    end
    seeds = [101,202,303];
    filesExist = true;
    for k = 1:numel(seeds)
        fileName = sprintf('seed_%d_compact.mat',seeds(k));
        if exist(fullfile(od,fileName),'file') == 0
            filesExist = false;
        end
    end
    if filesExist
        restore_compact;
        fprintf('Compact results restored. Original timing preserved.\n');
        return;
    end
    grid=linspace(1e-4,.8,801);
    w=ones(801,1)/800;
    w([1,end])=w([1,end])/2;
    if ~exist(fullfile(od,'protocol_data.mat'),'file')
        rng(2026101601,'twister');
        targets=struct([]);
        spec=[.07,3;.07,6;.13,3;.13,6;-.07,3;-.13,6];
        for j=1:6
            cc=c;
            cc.amplitude=spec(j,1);
            cc.steepness=spec(j,2);
            targets(j).id=sprintf('analytic_%02d',j);
            targets(j).kind='analytic';
            targets(j).F=target_curve(grid,cc);
            targets(j).amplitude=spec(j,1);
            targets(j).steepness=spec(j,2);
        end
        truth=struct([]);
        attempts=0;
        targettimer=tic;
        while numel(truth)<4
            x=.6*(2*rand(1,4)-1);
            [~,~,ok]=feed_decode(x,c);
            if ~ok
                continue;
            end
            cf=c;
            cf.RelTol=2e-9;
            cf.AbsTol=1e-13;
            cf.output_dt=5;
            cf.max_step=30;
            r=observe(x,cf,grid);
            attempts=attempts+1;
            if ~r.feasible
                continue;
            end
            k=numel(truth)+1;
            truth(k).x=x;
            truth(k).r=r;
            targets(6+k).id=sprintf('reachable_%02d',k);
            targets(6+k).kind='reachable';
            targets(6+k).F=r.F;
            targets(6+k).amplitude=NaN;
            targets(6+k).steepness=NaN;
        end
        target_seconds=toc(targettimer);
        save(fullfile(od,'hidden_target_truth.mat'),'truth','attempts','target_seconds');
        rng(2026101602,'twister');
        Xt=sample_controls(96,c);
        test=make_data(Xt,c,grid);
        save(fullfile(od,'protocol_data.mat'),'targets','test','c','grid','w','attempts','target_seconds');
        fprintf('PREPARED: test=96, target calls=%d, seconds=%.2f\n',attempts,toc(timer));
    else
        load(fullfile(od,'protocol_data.mat'),'targets','test','c','grid','w');
    end
    seeds=[101,202,303];
    for si=1:3
        maxNumCompThreads(1);
        seed=seeds(si);
        file=fullfile(od,sprintf('seed_%d.mat',seed));
        if ~exist(file,'file')
            res=run_seed(seed,c,grid,w,targets,test);
            save_seed(file,res);
        end
    end
    elapsed=toc(timer);
    save(fullfile(od,'run_elapsed.mat'),'elapsed');
    fprintf('PILOT COMPLETE: %.2f seconds wall time\n',elapsed);
end
function save_seed(path,res)
    save(path,'res','-v7.3');
end
function res=run_seed(seed,c,grid,w,targets,test)
    seedclock=tic;
    rng(seed,'twister');
    X=[zeros(1,4);sample_controls(159,c)];
    D=make_data(X,c,grid);
    assert(min(pdist2(D.X,test.X),[],'all')>1e-7,'Train/test overlap');
    t=tic;
    shared=fit_models(D,[],true);
    shared_fit_seconds=toc(t);
    t=tic;
    [Fpred,Mpred]=predict_models(shared,test.X,true);
    test_prediction_seconds=toc(t);
    valid=all(isfinite(test.F),2);
    e=Fpred(valid,:)-test.F(valid,:);
    tr=struct;
    tr.trajectory_rmse=sqrt((e.^2)*w);
    tr.max_error=max(abs(e),[],2);
    tr.test_indices=find(valid);
    tr.Fpred=Fpred;
    tr.Mpred=Mpred;
    tr.actual_feasible=test.feasible;
    tr.predicted_feasible=all(Mpred>=.8,2);
    tr.pca_rank=size(shared.B,2);
    tr.explained=shared.explained;
    tr.train_reconstruction_rmse=shared.reconstruction;
    tr.Xmol_mae=mean(abs(Mpred(:,1)-test.M(:,1)));
    tr.Xmass_mae=mean(abs(Mpred(:,2)-test.M(:,2)));
    tr.false_feasible=sum(tr.predicted_feasible & ~test.feasible);
    tr.false_infeasible=sum(~tr.predicted_feasible & test.feasible);
    fprintf('SEED %d: training %.2fs, feasible %d/160, rank %d, test trajectory RMSE %.5f\n',seed,sum(D.seconds),sum(D.feasible),tr.pca_rank,mean(tr.trajectory_rmse));
    runs=struct([]);
    idx=0;
    methods={'trajectory_gp','scalar_gp','sqp','random'};
    for ti=1:numel(targets)
        target=targets(ti);
        [base,ib]=best_score(D,target.F,w);
        baseline=struct('x',D.X(ib,:),'RMSE',base,'record',get_record(D,ib));
        for mi=1:numel(methods)
            name=methods{mi};
            t=tic;
            DD=D;
            fit_seconds=0;
            prediction_seconds=0;
            extra=[];
            if mi<=2
                for round=1:3
                    tfit=tic;
                    if mi==1 && round==1
                        model=shared;
                    else
                        model=fit_models(DD,target.F,mi==1);
                    end
                    fit_seconds=fit_seconds+toc(tfit);
                    [~,bestidx]=best_score(DD,target.F,w);
                    rng(seed+10000*ti+100*round,'twister');
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
                    [pred,pm]=predict_models(model,pool,mi==1);
                    if mi==1
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
                    extra=[extra;pool(chosen,:)];
                    
                end
            elseif mi==3
                rng(seed+10000*ti+777,'twister');
                DD=run_sqp(DD,target.F,w,c,grid,30);
                extra=DD.X(161:end,:);
            else
                rng(seed+10000*ti+888,'twister');
                new=make_data(sample_controls(30,c),c,grid);
                DD=append_data(DD,new);
                extra=new.X;
            end
            total_seconds=toc(t);
            assert(size(extra,1)==30);
            [score,bestidx]=best_score(DD,target.F,w);
            values=score_data(DD,target.F,w);
            hist=cummin(values);
            idx=idx+1;
            runs(idx).target=target;
            runs(idx).method=name;
            runs(idx).seed=seed;
            runs(idx).baseline=baseline;
            runs(idx).D=DD;
            runs(idx).best_x=DD.X(bestidx,:);
            runs(idx).best_record=get_record(DD,bestidx);
            runs(idx).RMSE=score;
            runs(idx).best10=hist(170);
            runs(idx).best20=hist(180);
            runs(idx).best30=hist(190);
            runs(idx).fit_seconds=fit_seconds;
            runs(idx).prediction_seconds=prediction_seconds;
            runs(idx).total_seconds=total_seconds;
            runs(idx).extra_ode_seconds=sum(DD.seconds(161:end));
            fprintf('SEED %d TARGET %s METHOD %s: base %.5f -> %.5f; %.2fs\n',seed,target.id,name,base,score,total_seconds);
        end
    end
    res=struct('seed',seed,'training',D,'test_metrics',tr,'shared_fit_seconds',shared_fit_seconds,'test_prediction_seconds',test_prediction_seconds,'runs',runs,'seconds',toc(seedclock));
end
function D=run_sqp(D,target,w,c,grid,budget)
    [~,ix]=sort(score_data(D,target,w));
    starts=D.X(ix(1:2),:);
    count=0;
    cap=0;
    for j=1:2
        cap=j*15;
        opt=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',15,'FiniteDifferenceStepSize',2e-3,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'MaxIterations',30);
        try
            fmincon(@obj,starts(j,:),[],[],[],[],-ones(1,4),ones(1,4),@cons,opt);
        catch ex
            if ~strcmp(ex.identifier,'pilot:budget')
                rethrow(ex);
            end
        end
    end
    while count<budget
        obj(sample_controls(1,c));
    end
    function y=obj(x)
        if count>=cap
            error('pilot:budget','budget');
        end
        new=make_data(x(:).',c,grid);
        D=append_data(D,new);
        count=count+1;
        y=score_data(new,target,w);
        if ~isfinite(y)
            y=10;
        end
    end
    function [g,h]=cons(x)
        [~,q]=feed_decode(x,c);
        g=[reshape(q./c.max_pump-1,[],1);sum(q,1).'/c.max_total_flow-1];
        h=[];
    end
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
