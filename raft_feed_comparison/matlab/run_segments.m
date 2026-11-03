function run_segments
% 使用相同初始方案和评价次数，比较不同进料分段。
    root = fileparts(fileparts(mfilename('fullpath')));
    base = fileparts(root);
    restoredefaultpath;
    addpath(fullfile(root,'matlab','core'),fullfile(root,'matlab'));
    maxNumCompThreads(1);
    od = fullfile(root,'outputs');
    diary(fullfile(od,'segments_log.txt'));
    for target = ["A16_K8","A16_K2"]
        for seed = [11 22 33]
            for variant = ["four","ten"]
                tag0 = sprintf('%s_%s_s%d',target,variant,seed);
                h = load(fullfile(base,'raft_high_budget','outputs',[tag0 '.mat']));
                c = h.run.c;
                c.budget = 450;
                c.sqp_start_budget = 90;
                for method = ["SQP","GP"]
                    tag = sprintf('segments_%s_%s',tag0,method);
                    dest = fullfile(od,[tag '.mat']);
                    if isfile(dest)
                        continue
                    end
                    fprintf('START %s %s\n',tag,string(datetime('now')));
                    timer = tic;
                    if method=="SQP"
                        [run,exceptions] = optimize(c,h.run.edges,h.run.initial,seed,tag,od);
                    else
                        [run,exceptions] = gp_search(c,h.run.edges,h.run.initial,seed,tag,od);
                    end
                    run.seconds = toc(timer);
                    run.target = target;
                    run.variant = variant;
                    run.method = method;
                    run.seed = seed;
                    save(dest,'run','exceptions','-v7');
                    fprintf('DONE %s RMSE %.9f calls %d seconds %.1f\n',tag,run.best_score,run.count,run.seconds);
                end
            end
        end
    end
    diary off
end
function x=encode(w)
    x=reshape(w(:,2:end).',1,[]);
end
function w=decode(x,K)
    t=reshape(x,K-1,2).';
    w=[1-sum(t,2),t];
end
function [A,b]=linear_constraints(c,edges)
    K=numel(edges)-1;
    D=2*(K-1);
    z=zeros(1,D);
    v0=reshape(decode(z,K).',[],1);
    M=zeros(2*K,D);
    for j=1:D
        e=z;
        e(j)=1;
        M(:,j)=reshape(decode(e,K).',[],1)-v0;
    end
    L=-eye(2*K);
    rhs=zeros(2*K,1);
    vm=c.p.MW./(1000*c.p.rho_monomer_g_mL);
    fac=reshape((c.dose.*vm./diff(edges)).',[],1);
    caps=repelem(c.max_pump,K);
    L=[L;diag(fac./caps)];
    rhs=[rhs;ones(2*K,1)];
    for k=1:K
        a=zeros(1,2*K);
        a(k)=fac(k)/c.max_total_flow;
        a(K+k)=fac(K+k)/c.max_total_flow;
        L=[L;a];
        rhs=[rhs;1];
        
    end
    % Preserve the exact early-stage logit range of the old controls.
    for species=0:1
        for j=1:2
            a=zeros(1,2*K);
            a(species*K+j)=1;
            a(species*K+3)=-exp(c.logit_bound);
            L=[L;a];
            rhs=[rhs;0];
            a=-a;
            a(species*K+3)=exp(-c.logit_bound);
            L=[L;a];
            rhs=[rhs;0];
            
        end
    end
    A=L*M;
    b=rhs-L*v0;
end
function [v,r,o]=score(x,c,edges)
    o=[];
    K=numel(edges)-1;
    w=decode(x,K);
    r=struct('score',NaN,'feasible',false,'RMSE',NaN,'max_error',NaN,'target_met',false,'Xmol',NaN,'Xmass',NaN,'Vmax',NaN,'peak_pump_mL_s',NaN,'message','','ode_attempted',false);
    if min(w,[],'all') < -1e-10
        v=2+sum(max(-w,0),'all');
        r.message='Negative dose fraction';
        r.score=v;
        return
    end
    w=max(w,0);
    w=w./sum(w,2);
    rates=c.dose.*w./diff(edges);
    q=rates.*(c.p.MW./(1000*c.p.rho_monomer_g_mL));
    r.peak_pump_mL_s=max(q,[],'all')*1000;
    if any(q>c.max_pump+1e-12,'all') || any(sum(q,1)>c.max_total_flow+1e-12)
        v=2+sum(max(q./c.max_pump-1,0),'all')+sum(max(sum(q,1)/c.max_total_flow-1,0));
        r.message='Pump cap exceeded';
        r.score=v;
        return
    end
    try
        cc=c;
        cc.feed_edges=edges;
        cc.feed_time=edges(end);
        r.ode_attempted=true;
        o=resolution_simulate(cc,rates);
        r.Xmol=o.Xmol(end);
        r.Xmass=o.Xmass(end);
        r.Vmax=max(o.V);
        violation=max(0,c.min_Xmol-r.Xmol)+max(0,c.min_Xmass-r.Xmass)+max(0,r.Vmax/c.max_volume-1);
        if violation>1e-8
            v=1+10*violation;
            r.message='Production or volume constraint';
        else
            [r.RMSE,r.max_error]=metrics(o,c,801);
            r.feasible=true;
            r.target_met=r.RMSE<=.02 && r.max_error<=.05;
            v=r.RMSE;
        end
    catch ex
        v=10;
        r.message=getReport(ex,'basic','hyperlinks','off');
    end
    r.score=v;
end
function [run,exceptions]=optimize(c,edges,initial,seed,tag,od)
    rng(seed,'twister');
    D=size(initial,2);
    K=numel(edges)-1;
    [A,b]=linear_constraints(c,edges);
    assert(max(A*initial.'-b,[],'all')<1e-9,'Initial pump/ratio constraints');
    X=zeros(c.budget,D);
    Y=nan(c.budget,1);
    records=[];
    count=0;
    odes=0;
    exceptions={};
    starts={};
    cap=c.budget;
    for j=1:3
        objective(initial(j,:));
    end
    start=0;
    while count<c.budget
        start=start+1;
        if start<=3
            x0=initial(start,:);
        else
            [~,ib]=min(Y(1:count));
            best=decode(X(ib,:),K);
            % A convex mixture with a shared legal coarse policy stays feasible.
            coarse=decode(initial(randi(3),:),K);
            alpha=.05+.25*rand;
            x0=encode((1-alpha)*best+alpha*coarse);
        end
        before=count;
        cap=min(c.budget,count+90);
        opts=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',90,...
        'FiniteDifferenceStepSize',2e-3,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'MaxIterations',200);
        exitflag=NaN;
        msg='';
        try
            [~,~,exitflag,info]=fmincon(@objective,x0,A,b,[],[],zeros(1,D),ones(1,D),[],opts);
            msg=info.message;
        catch ex
            if strcmp(ex.identifier,'resolution:budget')
                msg='Allocated evaluation cap reached';
            else
                msg=getReport(ex,'extended','hyperlinks','off');
                exceptions{end+1}=msg;
            end 
        end
        starts{end+1}=struct('index',start,'first_call',before+1,'last_call',count,'exitflag',exitflag,'message',msg);
        
        if count==before
            error('No objective evaluation in SQP start');
        end
    end
    [value,ix]=min(Y(1:count));
    run=struct('c',c,'edges',edges,'initial',initial,'X',X(1:count,:),'Y',Y(1:count),'records',records,'count',count,'ode_attempts',odes,'starts',{starts},'best_x',X(ix,:),'best_score',value,'best_record',records(ix),'best450',min(Y(1:min(450,count))));
    T=struct2table(records);
    T.evaluation=(1:count).';
    T.best_so_far=cummin(Y(1:count));
    writetable(T,fullfile(od,[tag '_history.csv']));
    function value=objective(x)
        if count>=cap
            error('resolution:budget','Evaluation cap reached');
        end
        [value,r]=score(x,c,edges);
        count=count+1;
        odes=odes+r.ode_attempted;
        X(count,:)=x(:).';
        Y(count)=value;
        if isempty(records)
            records=r;
        else
            records(count,1)=r;
        end
        if mod(count,150)==0
            fprintf('%s call %d best %.9f\n',tag,count,min(Y(1:count)));
            partial=struct('count',count,'X',X(1:count,:),'Y',Y(1:count),'records',records,'odes',odes);
            save(fullfile(od,[tag '_partial.mat']),'partial');
        end
    end
end
function [rmse,peak,grid,f,residual]=metrics(o,c,N)
    valid=o.Xmol>1e-8;
    [u,ix]=unique(o.Xmol(valid));
    F=o.Finstant(valid);
    F=F(ix);
    grid=linspace(1e-4,c.score_Xmol,N);
    assert(u(end)>=grid(end)-1e-8,'Insufficient production');
    f=interp1(u,F,grid,'linear','extrap');
    residual=f-target_curve(grid,c);
    rmse=sqrt(trapz(grid,residual.^2)/diff(grid([1 end])));
    peak=max(abs(residual));
end
function [run,exceptions]=gp_search(c,edges,initial,seed,tag,od)
    % 与SQP共用剂量份额坐标、线性可行域和初始三个物理方案。
    rng(seed,'twister');
    D=size(initial,2);
    K=numel(edges)-1;
    [A,b]=linear_constraints(c,edges);
    A=[A;-eye(D);eye(D)];
    b=[b;zeros(D,1);ones(D,1)];
    X=zeros(c.budget,D);
    Y=nan(c.budget,1);
    records=[];
    count=0;
    odes=0;
    exceptions={};
    chain=encode(ones(2,K)/K);
    assert(max(A*chain.'-b)<1e-9,'Interior sampling start is infeasible');
    for j=1:size(initial,1)
        add_point(initial(j,:));
    end
    while count<c.budget
        feasible=find([records.feasible]);
        [xu,ix]=unique(X(feasible,:),'rows','stable');
        yu=log(max(Y(feasible(ix)),1e-6));
        gp=fitrgp(xu,yu,'KernelFunction','ardmatern52','Standardize',true,'Sigma',1e-3,'ConstantSigma',true,'SigmaLowerBound',1e-6,'FitMethod','exact','PredictMethod','exact');
        [~,best]=min(Y(1:count));
        [globalPoints,chain]=sample_polytope(chain,A,b,1536);
        localPoints=zeros(512,D);
        for j=1:512
            alpha=.15*rand;
            localPoints(j,:)=(1-alpha)*X(best,:)+alpha*globalPoints(j,:);
        end
        candidates=[globalPoints;localPoints];
        [mu,sd]=predict(gp,candidates);
        sd=max(sd,1e-10);
        improvement=min(yu)-mu-.001;
        z=improvement./sd;
        ei=improvement.*normcdf(z)+sd.*normpdf(z);
        margins=min([[records.Xmol].'-c.min_Xmol,[records.Xmass].'-c.min_Xmass,1-[records.Vmax].'/c.max_volume],[],2);
        known=find(isfinite(margins));
        [cx,ci]=unique(X(known,:),'rows','stable');
        cg=fitrgp(cx,margins(known(ci)),'KernelFunction','ardmatern52','Standardize',true,'Sigma',1e-4,'ConstantSigma',true,'SigmaLowerBound',1e-7,'FitMethod','exact','PredictMethod','exact');
        [cm,cs]=predict(cg,candidates);
        ei=ei.*normcdf(cm./max(cs,1e-10));
        [~,ib]=max(ei);
        add_point(candidates(ib,:));
    end
    [value,ix]=min(Y);
    run=struct('c',c,'edges',edges,'initial',initial,'X',X,'Y',Y,'records',records,'count',count,'ode_attempts',odes,'best_x',X(ix,:),'best_score',value,'best_record',records(ix));
    function add_point(x)
        [value,r]=score(x,c,edges);
        count=count+1;
        X(count,:)=x;
        Y(count)=value;
        odes=odes+r.ode_attempted;
        if isempty(records)
            records=r;
        else
            records(count,1)=r;
        end
        if ~isempty(r.message) && value==10
            exceptions{end+1}=r.message;
        end
        if mod(count,50)==0
            fprintf('%s call %d best %.9f\n',tag,count,min(Y(1:count)));
            partial=struct('X',X(1:count,:),'Y',Y(1:count),'records',records,'odes',odes);
            save(fullfile(od,[tag '_partial.mat']),'partial');
        end
    end
end
function [points,x]=sample_polytope(x,A,b,n)
    % 沿随机方向求允许步长区间，在可行多面体内部移动。
    points=zeros(n,numel(x));
    for j=1:n
        for repeat=1:2
            direction=randn(numel(x),1);
            direction=direction/norm(direction);
            slope=A*direction;
            slack=b-A*x.';
            upper=min(slack(slope>1e-12)./slope(slope>1e-12));
            lower=max(slack(slope< -1e-12)./slope(slope< -1e-12));
            assert(upper>=lower-1e-9);
            step=lower+rand*(upper-lower);
            x=x+step*direction.';
        end
        points(j,:)=x;
    end
    assert(max(A*points.'-b,[],'all')<1e-8);
end
