function analyze_extension
% 按独立目标汇总误差，再查看三个重复之间的波动。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    fd=fullfile(root,'figures');
    old=fullfile(fileparts(root),'raft_trajectory_pilot','outputs');
    maxNumCompThreads(1);
    load(fullfile(od,'protocol_data.mat'),'targets','c','grid','w');
    load(fullfile(od,'hidden_target_truth.mat'),'truth','attempts','odecalls','target_seconds','logrows');
    load(fullfile(od,'run_elapsed.mat'),'elapsed','shared_cost');
    cf=c;
    cf.RelTol=2e-9;
    cf.AbsTol=1e-13;
    cf.output_dt=5;
    cf.max_step=30;
    dense=linspace(1e-4,.8,5001);
    rows=struct([]);
    details=struct([]);
    cache=struct([]);
    history=struct([]);
    seeds=[101 202 303];
    methods={'scalar_gp','trajectory_gp'};
    alltrain=[];
    allquery=[];
    match=true;
    querycalls=0;
    odecount=0;
    maxbalance=0;
    fine_secs=0;
    for si=1:3
        z=load(fullfile(old,sprintf('seed_%d_compact.mat',seeds(si))),'res');
        D=z.res.training;
        alltrain=[alltrain;D.X];
        for ti=1:8
            target=targets(ti);
            for mi=1:2
                z=load(fullfile(od,sprintf('run_%s_%d_%s.mat',target.id,seeds(si),methods{mi})),'run');
                run=z.run;
                match=match && isequal(D.X,run.initial_X) && isequaln(D.F,run.initial_F) && isequaln(D.M,run.initial_M);
                assert(size(run.extra.X,1)==30);
                allquery=[allquery;run.extra.X];
                querycalls=querycalls+30;
                odecount=odecount+sum(run.extra.pump_ok);
                rr=replay(run.best_x,target);
                rr.seed=seeds(si);
                rr.target=target.id;
                rr.method=methods{mi};
                rr.normal_RMSE=run.RMSE;
                rr.online_seconds=run.total_seconds;
                rr.fit_seconds=run.fit_seconds;
                rr.prediction_seconds=run.prediction_seconds;
                rr.ODE_seconds=run.extra_ode_seconds;
                rr.best10=run.best10;
                rr.best20=run.best20;
                rr.best30=run.best30;
                rr.additional_calls=30;
                rows=[rows,rr];
                o=cache(end).o;
                % 直接读取计数记录，包括使用缓存的情况。
                [o,~]=cached(run.best_x);
                f=profile(o,dense);
                details=[details,struct('target',target.id,'seed',seeds(si),'method',methods{mi},'x',run.best_x,'F',f,'t',o.t,'xi',o.Xmol,'Ftime',o.Finstant)];
                for ei=1:30
                    E=run.extra;
                    val=sqrt(((E.F(ei,:)-target.F).^2)*w);
                    history=[history,struct('target',target.id,'seed',seeds(si),'method',methods{mi},'additional_evaluation',ei,'RMSE',val,'feasible',E.feasible(ei),'pump_ok',E.pump_ok(ei),'seconds',E.seconds(ei),'x1',E.X(ei,1),'x2',E.X(ei,2),'x3',E.X(ei,3),'x4',E.X(ei,4))];
                end
            end
            rr=replay(run.baseline.x,target);
            rr.seed=seeds(si);
            rr.target=target.id;
            rr.method='database';
            rr.normal_RMSE=run.baseline.RMSE;
            rr.online_seconds=0;
            rr.fit_seconds=0;
            rr.prediction_seconds=0;
            rr.ODE_seconds=0;
            rr.best10=rr.fine_RMSE;
            rr.best20=rr.fine_RMSE;
            rr.best30=rr.fine_RMSE;
            rr.additional_calls=0;
            rows=[rows,rr];
        end
    end
    R=struct2table(rows);
    R=sortrows(R,{'target','seed','method'});
    writetable(R,fullfile(od,'design_results.csv'));
    writetable(struct2table(history),fullfile(od,'search_history_1440.csv'));
    % 先按目标汇总；随机重复属于同一目标，不视为新独立目标。
    targetrows=struct([]);
    names={'database','scalar_gp','trajectory_gp'};
    for ti=1:8
        for mi=1:3
            a=R(strcmp(R.target,targets(ti).id)&strcmp(R.method,names{mi}),:);
            targetrows=[targetrows,struct('target',targets(ti).id,'method',names{mi},'mean_RMSE',mean(a.fine_RMSE),'min_RMSE',min(a.fine_RMSE),'max_RMSE',max(a.fine_RMSE),'mean_max_error',mean(a.fine_max_error),'worst_max_error',max(a.fine_max_error),'met',sum(a.target_met),'runs',3,'mean_seconds',mean(a.online_seconds),'min_seconds',min(a.online_seconds),'max_seconds',max(a.online_seconds))];
        end
    end
    T=struct2table(targetrows);
    writetable(T,fullfile(od,'per_target_summary.csv'));
    oldR=readtable(fullfile(old,'design_results.csv'),'TextType','string');
    summary=struct([]);
    for group={'new8','old4','combined12','old_analytic6'}
        for mi=1:3
            name=names{mi};
            anew=R(strcmp(R.method,name),:);
            ao=oldR(oldR.method==string(name)&oldR.kind=="reachable",:);
            if strcmp(group{1},'old_analytic6')
                ao=oldR(oldR.method==string(name)&oldR.kind=="analytic",:);
            end
            switch group{1}
            case 'new8'
                vals=anew.fine_RMSE;
                peaks=anew.fine_max_error;
                met=anew.target_met;
                times=anew.online_seconds;
                n=8;
                ids=string(anew.target);
            case 'combined12'
                vals=[ao.fine_RMSE;anew.fine_RMSE];
                peaks=[ao.fine_max_error;anew.fine_max_error];
                met=[ao.target_met;anew.target_met];
                times=[ao.online_seconds;anew.online_seconds];
                n=12;
                ids=[ao.target;string(anew.target)];
            otherwise
                vals=ao.fine_RMSE;
                peaks=ao.fine_max_error;
                met=ao.target_met;
                times=ao.online_seconds;
                n=numel(unique(ao.target));
                ids=ao.target;
            end
            targetall=0;
            for id=unique(ids).'
                targetall=targetall+all(met(ids==id));
            end
            summary=[summary,struct('group',group{1},'method',name,'independent_targets',n,'runs',numel(vals),'mean_RMSE',mean(vals),'mean_max_error',mean(peaks),'worst_max_error',max(peaks),'met_runs',sum(met),'all_seeds_met_targets',targetall,'mean_online_seconds',mean(times))];
        end
    end
    S=struct2table(summary);
    writetable(S,fullfile(od,'group_summary.csv'));
    a=T(strcmp(T.method,'scalar_gp'),:);
    b=T(strcmp(T.method,'trajectory_gp'),:);
    assert(isequal(a.target,b.target));
    delta=a.mean_RMSE-b.mean_RMSE;
    rng(2026101903,'twister');
    boot=mean(delta(randi(8,8,10000)),1);
    paired=struct('independent_targets',8,'mean_absolute_improvement',mean(delta),'relative_improvement_ratio_of_means',mean(delta)/mean(a.mean_RMSE),'target_level_bootstrap95',prctile(boot,[2.5 97.5]),'target_mean_wins_1e4',sum(delta>1e-4),'target_mean_ties_1e4',sum(abs(delta)<=1e-4),'target_mean_losses_1e4',sum(delta< -1e-4));
    truthX=vertcat(truth.x);
    z=load(fullfile(old,'protocol_data.mat'),'test');
    oldtruth=load(fullfile(old,'hidden_target_truth.mat'),'truth');
    audit=struct('shared_initial_records_exact',match,'minimum_target_training_distance',min(pdist2(truthX,alltrain),[],'all'),'minimum_target_test_distance',min(pdist2(truthX,z.test.X),[],'all'),'minimum_target_old_target_distance',min(pdist2(truthX,vertcat(oldtruth.truth.x)),[],'all'),'minimum_target_new_search_distance',min(pdist2(truthX,allquery),[],'all'),'search_calls',querycalls,'search_actual_ODE',odecount,'all_final_feasible',all(R.feasible),'max_normal_fine_RMSE_difference',max(abs(R.normal_RMSE-R.fine_RMSE)),'maximum_balance_residual',max(R.balance),'fine_unique_calls',numel(cache),'fine_seconds',fine_secs);
    assert(match && querycalls==1440 && audit.minimum_target_training_distance>1e-6 && audit.minimum_target_test_distance>1e-6);
    cost=struct('new_training_calls',0,'reused_training_calls',480,'new_generation_proposals',attempts,'new_generation_ODE',odecalls,'generation_seconds',target_seconds,'search_calls',querycalls,'search_actual_ODE',odecount,'fine_unique_calls',numel(cache),'successful_production_ODE',odecalls+odecount+numel(cache),'failed_first_start_generation_ODE',1,'total_including_failed_first_start_ODE',1+odecalls+odecount+numel(cache),'fine_seconds',fine_secs,'successful_search_wall_seconds',elapsed,'online_sum_seconds',sum(R.online_seconds),'shared_initial_fit_seconds',sum([shared_cost.shared_fit_seconds]));
    report=struct('summary',summary,'paired_new8',paired,'audit',audit,'cost',cost,'shared_fit_cost',shared_cost);
    fid=fopen(fullfile(od,'results_summary.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
    fclose(fid);
    save(fullfile(od,'analysis.mat'),'R','T','S','paired','audit','cost','details','dense','targets','grid','cache','-v7.3');
    makeplots;
    disp(S);
    disp(paired);
    disp(audit);
    disp(cost);
    function rr=replay(x,target)
        [o,secs]=cached(x);
        fine_secs=fine_secs+secs;
        f=profile(o,dense);
        tar=interp1(grid,target.F,dense,'pchip');
        err=sqrt(trapz(dense,(f-tar).^2)/(dense(end)-dense(1)));
        peak=max(abs(f-tar));
        feas=o.Xmol(end)>=.8-1e-8 && o.Xmass(end)>=.8-1e-8 && max(o.V)<=c.max_volume+1e-8;
        rr=struct('fine_RMSE',err,'fine_max_error',peak,'feasible',feas,'target_met',feas && err<=.02 && peak<=.05,'Xmol',o.Xmol(end),'Xmass',o.Xmass(end),'Vmax',max(o.V),'balance',max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]));
    end
    function [o,secs]=cached(x)
        secs=0;
        for k=1:numel(cache)
            if max(abs(cache(k).x-x))<1e-10
                o=cache(k).o;
                return;
            end
        end
        t=tic;
        [rates,~,ok,~,cc]=feed_decode(x,cf);
        assert(ok);
        o=simulate_reactor(cc,rates);
        secs=toc(t);
        cache=[cache,struct('x',x,'o',o)];
    end
    function makeplots
        set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',11);
        colors=[.16 .42 .68;.85 .35 .12];
        f=figure('Visible','off','Color','w','Position',[80 80 1180 780]);
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        nexttile;
        hold on;
        for mm=1:2
            aa=T(strcmp(T.method,methods{mm}),:);
            xx=(5:12)+(mm-1.5)*.13;
            errorbar(xx,aa.mean_RMSE,aa.mean_RMSE-aa.min_RMSE,aa.max_RMSE-aa.mean_RMSE,'o-','Color',colors(mm,:),'LineWidth',1.2,'MarkerFaceColor',colors(mm,:));
        end
        yline(.02,'k--');
        xticks(5:12);
        xlabel('Independent hidden target');
        ylabel('Tracking RMSE');
        title('(a) Target means and seed ranges');
        legend('Scalar GP','Trajectory GP','RMSE threshold','Location','best');
        set(gca,'XGrid','on','YGrid','on')
        nexttile;
        hold on;
        for mm=1:2
            aa=sortrows(R(strcmp(R.method,methods{mm}),:),{'target','seed'});
            scatter(repelem((5:12)',3)+(mm-1.5)*.15,aa.fine_max_error,24,colors(mm,:),'filled');
        end
        yline(.05,'k--');
        xticks(5:12);
        xlabel('Independent hidden target');
        ylabel('Maximum absolute error');
        title('(b) All three seeds retained');
        set(gca,'XGrid','on','YGrid','on')
        nexttile;
        bar(5:12,100*delta./a.mean_RMSE,'FaceColor',colors(2,:));
        yline(0,'k-');
        xticks(5:12);
        xlabel('Independent hidden target');
        ylabel('RMSE reduction relative to scalar GP (%)');
        title('(c) Per-target improvement');
        set(gca,'XGrid','on','YGrid','on')
        nexttile;
        costs=zeros(2,3);
        for mm=1:2
            aa=R(strcmp(R.method,methods{mm}),:);
            costs(mm,:)=[mean(aa.ODE_seconds),mean(aa.fit_seconds),mean(aa.online_seconds-aa.ODE_seconds-aa.fit_seconds)];
        end
        bar(costs,'stacked');
        xticklabels({'Scalar GP','Trajectory GP'});
        ylabel('Mean online time per run (s)');
        legend('Forward ODE','GP fitting','Candidate selection and other','Location','northwest');
        title('(d) New eight targets: equal 30-call budget');
        set(gca,'XGrid','on','YGrid','on')
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(f);
        exportgraphics(f,fullfile(fd,'extension_summary.png'),'Resolution',360);
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(f);
        exportgraphics(f,fullfile(fd,'extension_summary.pdf'),'ContentType','vector');
        savefig(f,fullfile(fd,'extension_summary.fig'));
        close(f);
        f=figure('Visible','off','Color','w','Position',[80 80 1400 700]);
        tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
        for ti=1:8
            nexttile;
            hold on;
            plot(grid,targets(ti).F,'k--','LineWidth',1.4);
            for mm=1:2
                ix=find(strcmp({details.target},targets(ti).id)&[details.seed]==202&strcmp({details.method},methods{mm}));
                plot(dense,details(ix).F,'Color',colors(mm,:),'LineWidth',1.1);
            end
            title(strrep(targets(ti).id,'_',' '));
            xlabel('Polymerization progress');
            ylabel('Instantaneous St fraction');
            set(gca,'XGrid','on','YGrid','on')
            if ti==1
                legend('Target','Scalar GP','Trajectory GP','Location','best');
            end
        end
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(f);
        exportgraphics(f,fullfile(fd,'extension_all_trajectories.png'),'Resolution',360);
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(f);
        exportgraphics(f,fullfile(fd,'extension_all_trajectories.pdf'),'ContentType','vector');
        savefig(f,fullfile(fd,'extension_all_trajectories.fig'));
        close(f);
    end
end
function f=profile(o,grid)
    valid=o.Xmol>1e-8;
    [u,ix]=unique(o.Xmol(valid));
    F=o.Finstant(valid);
    f=interp1(u,F(ix),grid,'linear');
end
