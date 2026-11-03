function run_resolution(mode)
% 计算较高预算下的分段进料结果，并保留每次搜索记录。
    % 四段与十段进料配对比较；动力学不变，不拟合GP。
    if nargin<1
        mode='all';
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    core=fullfile(root,'matlab','core');
    savedpath=path;
    cl=onCleanup(@()path(savedpath));
    
    restoredefaultpath;
    addpath(core,'-begin');
    addpath(fullfile(root,'matlab'),'-begin');
    assert(strcmp(which('reactor_rhs'),fullfile(core,'reactor_rhs.m')));
    maxNumCompThreads(1);
    od=fullfile(root,'outputs');
    diary(fullfile(od,['run_' mode '.log']));
    dc=onCleanup(@()diary('off'));
    
    fprintf('START %s MODE %s\n',char(datetime('now','TimeZone','Asia/Shanghai')),mode);
    if any(strcmp(mode,{'all','check'}))
        mapping_checks(root);
    end
    if any(strcmp(mode,{'all','search'}))
        for target=["A16_K8","A16_K2"]
            for seed=[11 22 33]
                h=load(fullfile(root,'historical',sprintf('intervention_%s_late_tail_s%d.mat',target,seed)));
                c=h.run.c;
                c.budget=1500;
                rng(seed,'twister');
                physical={h.run.best_x,neutral_control(c),sample_controls(1,c)};
                [r0,~,~,~,cc]=feed_decode(physical{1},c);
                %#ok<ASGLU>
                ed4=[0 2400 4800 7200 21600];
                W4=zeros(3,6);
                W10=zeros(3,18);
                for j=1:3
                    [rates,~,ok]=feed_decode(physical{j},c);
                    assert(ok);
                    w4=rates.*diff(ed4)./c.dose;
                    w10=[w4(:,1:3),repmat(w4(:,4)/6,1,6),zeros(2,1)];
                    W4(j,:)=encode(w4);
                    W10(j,:)=encode(w10);
                end
                for variant=["four","ten"]
                    if variant=="four"
                        edges=ed4;
                        initial=W4;
                    else
                        edges=0:2400:24000;
                        initial=W10;
                    end
                    tag=sprintf('%s_%s_s%d',target,variant,seed);
                    dest=fullfile(od,[tag '.mat']);
                    if exist(dest,'file')
                        fprintf('REUSE %s\n',tag);
                        continue;
                    end
                    fprintf('BEGIN %s\n',tag);
                    timer=tic;
                    [run,exceptions]=optimize(c,edges,initial,seed,tag,od);
                    run.seconds=toc(timer);
                    run.target=target;
                    run.variant=variant;
                    run.seed=seed;
                    run.initial_physical_controls=physical;
                    save(dest,'run','exceptions','-v7');
                    fprintf('COMPLETE %s calls %d ODE %d best %.9f time %.2f\n',tag,run.count,run.ode_attempts,run.best_score,run.seconds);
                end
            end
        end
    end
    if any(strcmp(mode,{'all','analyze'}))
        analyze(root);
    end
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
    % 保留早期控制原有的logit取值范围。
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
            % 与共同可行粗分段方案作凸组合，可保持此处的进料约束。
            coarse=decode(initial(randi(3),:),K);
            alpha=.05+.25*rand;
            x0=encode((1-alpha)*best+alpha*coarse);
        end
        before=count;
        cap=min(c.budget,count+500);
        opts=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',500,...
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
function mapping_checks(root)
    od=fullfile(root,'outputs');
    dest=fullfile(od,'mapping_checks.csv');
    if exist(dest,'file')
        return;
    end
    rows={};
    for target=["A16_K8","A16_K2"]
        for seed=[11 22 33]
            h=load(fullfile(root,'historical',sprintf('intervention_%s_late_tail_s%d.mat',target,seed)));
            c=h.run.c;
            [rates,~,~,~,cc]=feed_decode(h.run.best_x,c);
            old=simulate_reactor(cc,rates);
            new=resolution_simulate(cc,rates);
            assert(isequal(old.y,new.y),'Integration wrapper changed historical trajectory');
            w4=rates.*[2400 2400 2400 14400]./c.dose;
            w10=[w4(:,1:3),repmat(w4(:,4)/6,1,6),zeros(2,1)];
            [v4,r4]=score(encode(w4),c,[0 2400 4800 7200 21600]);
            [v10,r10]=score(encode(w10),c,0:2400:24000);
            assert(abs(v4-v10)<1e-4 && abs(v4-h.run.best_score)<1e-8);
            rows{end+1}=struct('target',target,'seed',seed,'historical_RMSE',h.run.best_score,'four_RMSE',v4,'mapped_ten_RMSE',v10,'mapping_difference',abs(v4-v10),'four_Xmol',r4.Xmol,'ten_Xmol',r10.Xmol,'ode_calls',4);
            
        end
    end
    writetable(struct2table(vertcat(rows{:})),dest);
    fprintf('MAPPING CHECKS passed: 24 diagnostic ODE calls\n');
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
function analyze(root)
    od=fullfile(root,'outputs');
    files=dir(fullfile(od,'A16_K*_s*.mat'));
    rows={};
    feedrows={};
    stored=struct;
    for fi=1:numel(files)
        if contains(files(fi).name,'partial')
            continue;
        end
        h=load(fullfile(od,files(fi).name));
        if ~isfield(h,'run')
            continue;
        end;
        r=h.run;
        tag=erase(files(fi).name,'.mat');
        dest=fullfile(od,[tag '_fine.mat']);
        if exist(dest,'file')
            q=load(dest);
            o=q.o;
            rr=q.rr;
            fine=q.fine;
        else
            fine=r.c;
            fine.RelTol=2e-9;
            fine.AbsTol=1e-13;
            fine.output_dt=5;
            fine.max_step=30;
            [~,rr,o]=score(r.best_x,fine,r.edges);
            save(dest,'o','rr','fine');
        end
        assert(~isempty(o),'Best policy fine integration failed');
        [rmse,peak,grid,f,residual]=metrics(o,r.c,5001);
        [~,mi]=max(abs(residual));
        w=decode(r.best_x,numel(r.edges)-1);
        rates=r.c.dose.*w./diff(r.edges);
        q=rates.*(r.c.p.MW./(1000*r.c.p.rho_monomer_g_mL))*1000;
        dose_error=max(abs(sum(rates.*diff(r.edges),2)-r.c.dose));
        balance=max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]);
        feasible=rr.feasible && min(w,[],'all')>=-1e-10 && dose_error<1e-9;
        rows{end+1}=struct('target',r.target,'configuration',r.variant,'seed',r.seed,'calls',r.count,'ode_attempts',r.ode_attempts,'seconds',r.seconds,'best450',r.best450,'search_RMSE',r.best_score,'fine_RMSE',rmse,'max_error',peak,'peak_progress',grid(mi),'peak_signed_error',residual(mi),'target_met',feasible && rmse<=.02 && peak<=.05,'feasible',feasible,'Xmol',rr.Xmol,'Xmass',rr.Xmass,'Vmax',rr.Vmax,'max_pump_mL_s',max(q,[],'all'),'max_total_mL_s',max(sum(q,1)),'fine_change',abs(rr.RMSE-r.best_score),'grid_change',abs(rmse-rr.RMSE),'balance',balance,'dose_error_mol',dose_error,'starts',numel(r.starts));
        
        for j=1:size(q,2)
            feedrows{end+1}=struct('target',r.target,'configuration',r.variant,'seed',r.seed,'start_s',r.edges(j),'end_s',r.edges(j+1),'St_mL_s',q(1,j),'GMA_mL_s',q(2,j));
        end 
        trajectory=table(grid.',f.',target_curve(grid,r.c).',residual.','VariableNames',{'progress','St_fraction','target','residual'});
        writetable(trajectory,fullfile(od,[tag '_trajectory.csv']));
        stored.(matlab.lang.makeValidName(tag))=struct('run',r,'o',o,'grid',grid,'f',f,'q',q);
    end
    R=struct2table(vertcat(rows{:}));
    assert(height(R)==12);
    writetable(R,fullfile(od,'results.csv'));
    writetable(struct2table(vertcat(feedrows{:})),fullfile(od,'feed_schedules.csv'));
    ar={};
    for target=["A16_K8","A16_K2"]
        for variant=["four","ten"]
            t=R(R.target==target & R.configuration==variant,:);
            ar{end+1}=struct('target',target,'configuration',variant,'runs',height(t),'mean_RMSE',mean(t.fine_RMSE),'std_RMSE',std(t.fine_RMSE),'best_RMSE',min(t.fine_RMSE),'mean_max_error',mean(t.max_error),'worst_max_error',max(t.max_error),'successes',sum(t.target_met),'mean_seconds',mean(t.seconds),'mean_best450',mean(t.best450));
            
        end
    end
    A=struct2table(vertcat(ar{:}));
    writetable(A,fullfile(od,'aggregate.csv'));
    save(fullfile(od,'analysis.mat'),'R','A','stored','-v7');
    S=struct('runs',height(R),'budget_calls',sum(R.calls),'search_ode_attempts',sum(R.ode_attempts),'mapping_ODE_calls',24,'fine_ODE_calls',12,'total_ODE_attempts',sum(R.ode_attempts)+36,'max_fine_change',max(R.fine_change),'max_grid_change',max(R.grid_change),'max_balance',max(R.balance),'max_dose_error',max(R.dose_error_mol),'all_feasible',all(R.feasible),'twenty_stage_runs',0,'GP_training_runs',0);
    fid=fopen(fullfile(od,'verification.json'),'w');
    fprintf(fid,'%s',jsonencode(S,PrettyPrint=true));
    fclose(fid);
    make_figures(R,stored,root);
    disp(A);
    disp(S);
end
function make_figures(R,S,root)
    f=figure('Visible','off','Color','w','Position',[100 100 1100 780]);
    tiledlayout(2,2,'TileSpacing','compact');
    ids=["A16_K8","A16_K2"];
    for ti=1:2
        target=ids(ti);
        nexttile(ti);
        hold on;
        plot([0 .8],[0 0],'Visible','off');
        for variant=["four","ten"]
            z=R(R.target==target & R.configuration==variant,:);
            [~,ib]=min(z.fine_RMSE);
            tag=sprintf('%s_%s_s%d',target,variant,z.seed(ib));
            s=S.(matlab.lang.makeValidName(tag));
            if variant=="four"
                plot(s.grid,target_curve(s.grid,s.run.c),'k--','DisplayName','Target','LineWidth',1.5);
            end
            plot(s.grid,s.f,'LineWidth',1.5,'DisplayName',char(variant));
        end
        ylim([.45 .95]);
        xlabel('Polymerization progress');
        ylabel('Instantaneous St fraction');
        title(strrep(target,'_',' / '));
        legend('Location','best');
        grid on;
        set(gca,'FontName','Times New Roman');
        nexttile(ti+2);
        hold on;
        for variant=["four","ten"]
            z=R(R.target==target & R.configuration==variant,:);
            curves=zeros(3,1500);
            for j=1:3
                tag=sprintf('%s_%s_s%d',target,variant,z.seed(j));
                s=S.(matlab.lang.makeValidName(tag));
                curves(j,:)=cummin(s.run.Y);
            end
            plot(1:1500,mean(curves,1),'LineWidth',1.5,'DisplayName',char(variant));
        end
        yline(.02,'k:','RMSE threshold','HandleVisibility','off');
        xlabel('New scoring calls');
        ylabel('Mean best score');
        legend('Location','best');
        grid on;
        set(gca,'FontName','Times New Roman');
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','resolution_comparison.png'),'Resolution',360);
    savefig(f,fullfile(root,'figures','resolution_comparison.fig'));
    close(f);
    f=figure('Visible','off','Color','w','Position',[100 100 1100 560]);
    tiledlayout(1,2);
    for ti=1:2
        z=R(R.target==ids(ti) & R.configuration=="ten",:);
        [~,ib]=min(z.fine_RMSE);
        tag=sprintf('%s_ten_s%d',ids(ti),z.seed(ib));
        s=S.(matlab.lang.makeValidName(tag));
        nexttile;
        hold on;
        stairs(s.run.edges/3600,[s.q(1,:),s.q(1,end)],'LineWidth',1.4,'DisplayName','St');
        stairs(s.run.edges/3600,[s.q(2,:),s.q(2,end)],'LineWidth',1.4,'DisplayName','GMA');
        xlabel('Time / h');
        ylabel('Feed rate / (mL s^{-1})');
        title(strrep(ids(ti),'_',' / '));
        legend;
        grid on;
        set(gca,'FontName','Times New Roman');
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','ten_stage_feeds.png'),'Resolution',360);
    savefig(f,fullfile(root,'figures','ten_stage_feeds.fig'));
    close(f);
end
