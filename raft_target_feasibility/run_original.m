function run_original(seed, target_indices)
% 从共同训练记录出发，在原三段进料范围内继续搜索。
    if nargin<2
        target_indices=1:6;
    end
    root=fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'core'));
    maxNumCompThreads(1);
    source=fullfile(fileparts(root),'raft_trajectory_pilot','outputs');
    p=load(fullfile(source,'protocol_data.mat'));
    s=load(fullfile(source,sprintf('seed_%d_compact.mat',seed)));
    D=s.res.training;
    c=p.c;
    grid=p.grid;
    w=p.w;
    od=fullfile(root,'outputs');
    if ~isfolder(od)
        mkdir(od);
    end
    for ti=target_indices
        dest=fullfile(od,sprintf('original_t%d_seed%d.mat',ti,seed));
        if isfile(dest)
            continue;
        end
        target=p.targets(ti);
        rng(seed+10000*ti+777,'twister');
        rv=sqrt(((D.F-target.F).^2)*w);
        pv=max(abs(D.F-target.F),[],2);
        rv(~D.feasible)=Inf;
        pv(~D.feasible)=Inf;
        [bestR,ir]=min(rv);
        [bestM,im]=min(max(rv/.02,pv/.05));
        xr=D.X(ir,:);
        xm=D.X(im,:);
        startsX=[xr;sample_controls(2,c)];
        count=0;
        odes=0;
        cap=0;
        cachex=[];
        obj=[Inf Inf];
        cv=[];
        history=nan(900,6);
        controls=zeros(900,4);
        errors={};
        exits={};
        timer=tic;
        fprintf('START t%d seed%d\n',ti,seed);
        for phase=1:2
            for st=1:3
                cap=(phase-1)*450+st*150;
                start=startsX(st,:);
                if phase==2&&st==1
                    start=xm;
                end
                options=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',150,'MaxIterations',100,'FiniteDifferenceStepSize',2e-3,'OptimalityTolerance',1e-6,'StepTolerance',1e-7);
                try
                    [~,~,flag,info]=fmincon(@objective,start,[],[],[],[],-ones(1,4),ones(1,4),@constraints,options);
                    msg=info.message;
                catch ex
                    flag=-99;
                    msg=ex.message;
                    if ~strcmp(ex.identifier,'diag:budget')
                        errors{end+1}=getReport(ex);
                    end
                end
                exits{end+1}=struct('phase',phase,'start',st,'flag',flag,'message',msg,'count',count);
                while count<cap
                    evaluate(sample_controls(1,c));
                end
            end
        end
        candidates=[xr;xm];
        finecell=cell(1,2);
        for k=1:2
            cf=c;
            cf.RelTol=5e-10;
            cf.AbsTol=2e-14;
            cf.output_dt=2.5;
            cf.max_step=15;
            [rates,q,ok,~,cc]=feed_decode(candidates(k,:),cf);
            o=simulate_reactor(cc,rates);
            gg=linspace(1e-4,.8,10001);
            valid=o.Xmol>1e-8;
            [u,ix]=unique(o.Xmol(valid));
            f=o.Finstant(valid);
            ff=interp1(u,f(ix),gg);
            truth=interp1(grid,target.F,gg,'pchip');
            e=ff-truth;
            rr=sqrt(trapz(gg,e.^2)/(.8-1e-4));
            [pk,ip]=max(abs(e));
            feas=ok&&o.Xmol(end)>=.8&&o.Xmass(end)>=.8&&max(o.V)<=c.max_volume;
            finecell{k}=struct('rmse',rr,'peak',pk,'peak_progress',gg(ip),'start_error',abs(e(1)),'end_error',abs(e(end)),'feasible',feas,'passed',feas&&rr<=.02&&pk<=.05,'Xmol',o.Xmol(end),'Xmass',o.Xmass(end),'maxV',max(o.V),'grid',gg,'F',ff,'target',truth,'q',q,'checks',o.checks);
        end
        fine=[finecell{:}];
        result=struct('target',target,'seed',seed,'candidates',candidates,'fine',fine,'count',count,'odes',odes,'history',history,'controls',controls,'errors',{errors},'exits',{exits},'seconds',toc(timer),'initial_best_rmse',min(rv));
        save(dest,'result','-v7');
        fprintf('DONE t%d seed%d RMSE %.6f peak %.6f PASS%d / RMSE %.6f peak %.6f PASS%d calls%d\n',ti,seed,fine(1).rmse,fine(1).peak,fine(1).passed,fine(2).rmse,fine(2).peak,fine(2).passed,count);
    end
    function evaluate(x)
        x=x(:).';
        if isequal(x,cachex)
            return;
        end
        if count>=cap
            error('diag:budget','budget');
        end
        count=count+1;
        cachex=x;
        controls(count,:)=x;
        obj=[10 1000];
        [rates,q,ok,~,cc]=feed_decode(x,c);
        cv=[reshape(q./c.max_pump-1,[],1);sum(q,1).'/c.max_total_flow-1];
        if ~ok
            history(count,:)=[NaN NaN 0 0 bestR bestM];
            return;
        end
        try
            odes=odes+1;
            o=simulate_reactor(cc,rates);
            feas=o.Xmol(end)>=.8&&o.Xmass(end)>=.8&&max(o.V)<=c.max_volume;
            if o.Xmol(end)>=.8
                valid=o.Xmol>1e-8;
                [u,ix]=unique(o.Xmol(valid));
                f=o.Finstant(valid);
                ff=interp1(u,f(ix),grid);
                err=ff-target.F;
                rr=sqrt((err.^2)*w);
                pk=max(abs(err));
                mm=max(rr/.02,pk/.05);
                if feas
                    obj=[rr mm];
                    if rr<bestR
                        bestR=rr;
                        xr=x;
                    end;
                    if mm<bestM
                        bestM=mm;
                        xm=x;
                    end
                end
                history(count,:)=[rr pk feas o.Xmol(end) bestR bestM];
            end
        catch ex
            errors{end+1}=ex.message;
        end
    end
    function v=objective(x)
        evaluate(x);
        v=obj(phase);
    end
    function [g,h]=constraints(x)
        [~,q]=feed_decode(x,c);
        g=[reshape(q./c.max_pump-1,[],1);sum(q,1).'/c.max_total_flow-1];
        h=[];
    end
end
