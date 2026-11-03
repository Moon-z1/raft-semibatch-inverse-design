function run_expanded(seed)
% 同时调整初始组成和全程进料分配，检查目标能否达标。
    % 独立的操作条件×评分预算对照，不读入旧试验的成功控制。
    root=fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'core'));
    maxNumCompThreads(1);
    od=fullfile(root,'outputs');
    if ~isfolder(od)
        mkdir(od);
    end
    ids=["A07_K3","A07_K6","A13_K3","A13_K6","D07_K3","D13_K6"];
    for id=ids
        for expanded=1:1
            dest=fullfile(od,sprintf('%s_op%d_seed%d.mat',id,expanded,seed));
            if isfile(dest)
                continue;
            end
            rng(seed,'twister');
            c=study_config;
            if id~="constant"
                v=sscanf(char(extractAfter(id,1)),'%d_K%d');
                c.amplitude=v(1)/100;
                c.steepness=v(2);
                if startsWith(id,"D")
                    c.amplitude=-c.amplitude;
                end;
            end
            c.target_id=char(id);
            total=c.n0+c.dose;
            vm=c.p.MW./(1000*c.p.rho_monomer_g_mL);
            if expanded
                K=10;
                edges=0:2400:24000;
            else
                K=3;
                edges=0:2400:7200;
            end
            N=K+1;
            D=2*N;
            dt=diff(edges);
            E=zeros(3,D);
            E(1,1:N)=1;
            E(2,N+1:end)=1;
            E(3,1)=total(1)*vm(1)/c.V0;
            E(3,N+1)=total(2)*vm(2)/c.V0;
            eq=ones(3,1);
            if ~expanded
                E(4,1)=1;
                eq(4)=c.n0(1)/total(1);
            end
            A=zeros(3*K,D);
            b=ones(3*K,1);
            for k=1:K
                a=total(1)*vm(1)/dt(k);
                g=total(2)*vm(2)/dt(k);
                A(k,1+k)=a/c.max_pump(1);
                A(K+k,N+1+k)=g/c.max_pump(2);
                A(2*K+k,1+k)=a/c.max_total_flow;
                A(2*K+k,N+1+k)=g/c.max_total_flow;
            end
            if ~expanded
                % 原三段softmax的logit_bound=4可用剂量比线性不等式精确表达。
                for sp=0:1
                    for k=1:K-1
                        row=zeros(1,D);
                        row(sp*N+k+1)=1;
                        row(sp*N+N)=-exp(c.logit_bound);
                        A(end+1,:)=row;
                        b(end+1)=0;
                        row=zeros(1,D);
                        row(sp*N+N)=1;
                        row(sp*N+k+1)=-exp(c.logit_bound);
                        A(end+1,:)=row;
                        b(end+1)=0;
                    end
                end
            end
            lb=zeros(D,1);
            ub=ones(D,1);
            w=zeros(2,N);
            w(:,1)=c.n0./total;
            w(:,2:4)=repmat(c.dose./total/3,1,3);
            x0=reshape(w.',[],1);
            % 两种操作配置从同一个物理进料方案出发，新增后段初值为零。
            directions=null(E);
            G=[A;-eye(D);eye(D)];
            h=[b;zeros(D,1);ones(D,1)];
            count=0;
            odes=0;
            best=Inf;
            bestx=x0;
            history=nan(450,8);
            controls=zeros(450,D);
            cachex=[];
            cacheobj=[];
            cachec=[];
            errors={};
            starts={};
            checkpoints={};
            cap=90;
            timer=tic;
            st=0;
            fprintf('START %s op%d seed%d\n',id,expanded,seed);
            while count<450
                st=st+1;
                cap=min(450,count+90);
                if count<450
                    cap=min(cap,450);
                end
                if st>1
                    x0=bestx;
                    if mod(st,3)==0
                        % 随机线性规划给出合法方向，混合得到可复现的新起点。
                        objectiveDirection=randn(D,1);
                        [vertex,~,flag]=linprog(objectiveDirection,A,b,E,eq,lb,ub,optimoptions('linprog','Display','none'));
                        if flag>0
                            x0=.8*x0+.2*vertex;
                        end
                    end
                end
                before=count;
                options=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',90,...
                'MaxIterations',150,'FiniteDifferenceStepSize',1e-4,'OptimalityTolerance',1e-6,...
                'ConstraintTolerance',1e-7,'StepTolerance',1e-8);
                try
                    [~,~,flag,info]=fmincon(@objective,x0,A,b,E,eq,lb,ub,@nonlin,options);
                    message=info.message;
                catch ex
                    flag=-99;
                    message=ex.message;
                    if ~strcmp(ex.identifier,'trial:budget')
                        errors{end+1}=getReport(ex,'extended','hyperlinks','off');
                    end
                end
                starts{end+1}=struct('first',before+1,'last',count,'flag',flag,'message',message);
                if count==450&&isempty(checkpoints)
                    checkpoints{1}=struct('budget',450,'x',bestx,'search_merit',best,'elapsed',toc(timer));
                end
                if count==before
                    error('No scoring progress');
                end
                if mod(st,5)==0
                    fprintf('%s op%d seed%d calls%d merit%.5f\n',id,expanded,seed,count,best);
                end
            end
            checkpoints{2}=struct('budget',450,'x',bestx,'search_merit',best,'elapsed',toc(timer));
            run=struct('id',id,'expanded',expanded,'seed',seed,'c',c,'edges',edges,'total',total,...
            'A',A,'b',b,'E',E,'eq',eq,'checkpoints',{checkpoints},'calls',count,'odes',odes,...
            'history',history,'controls',controls,'starts',{starts},'errors',{errors},'seconds',toc(timer));
            save(dest,'run','-v7');
            fprintf('DONE %s op%d seed%d time%.1f\n',id,expanded,seed,run.seconds);
        end
    end
    function cc=decode(x)
        amount=reshape(x,N,2).'.*total;
        cc=c;
        cc.n0=amount(:,1);
        cc.dose=total-cc.n0;
        cc.feed_edges=edges;
        cc.feed_time=edges(end);
        cc.test_rates=amount(:,2:end)./dt;
    end
    function evaluate(x)
        if isequal(x,cachex)
            return;
        end
        if count>=cap
            error('trial:budget','Local scoring allocation reached');
        end
        cachex=x;
        count=count+1;
        controls(count,:)=x.';
        cc=decode(x);
        cacheobj=1e4;
        cachec=ones(205,1)*10;
        try
            odes=odes+1;
            o=resolution_simulate(cc,cc.test_rates);
            [u,ix]=unique(o.Xmol);
            grid=linspace(1e-4,.8,801);
            err=interp1(u,o.Finstant(ix),grid,'linear','extrap')-target_curve(grid,c);
            rmse=sqrt(trapz(grid,err.^2)/(.8-1e-4));
            peak=max(abs(err));
            gg=linspace(1e-4,.8,101);
            ee=interp1(u,o.Finstant(ix),gg,'linear','extrap')-target_curve(gg,c);
            cacheobj=(rmse/.02)^2;
            cachec=[ee(:)/.05-1;-ee(:)/.05-1;(.80001-o.Xmol(end))/.8;(.80001-o.Xmass(end))/.8;max(o.V)/c.max_volume-1];
            linear=max([0;A*x-b;abs(E*x-eq);-x;x-1]);
            feasible=o.Xmol(end)>=.8&&o.Xmass(end)>=.8&&max(o.V)<=c.max_volume&&linear<1e-8;
            merit=max(rmse/.02,peak/.05)+100*double(~feasible);
            if merit<best
                best=merit;
                bestx=x;
            end
            history(count,:)=[rmse peak o.Xmol(end) o.Xmass(end) max(o.V) linear merit best];
        catch ex
            errors{end+1}=ex.message;
        end
    end
    function v=objective(x)
        evaluate(x);
        v=cacheobj;
    end
    function [g,h0]=nonlin(x)
        evaluate(x);
        g=cachec;
        h0=[];
    end
end
