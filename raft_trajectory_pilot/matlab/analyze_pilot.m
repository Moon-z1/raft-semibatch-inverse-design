function analyze_pilot
% 对最终进料进行精细积分，比较预测结果和实际跟踪误差。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    fd=fullfile(root,'figures');
    load(fullfile(od,'protocol_data.mat'),'targets','test','c','grid','w','attempts','target_seconds');
    seeds=[101,202,303];
    allruns=struct([]);
    testrows=struct([]);
    seedrows=struct([]);
    errors=cell(1,3);
    fine=struct([]);
    counter=0;
    replayseconds=0;
    finecache=struct([]);
    details=struct([]);
    cf=c;
    cf.RelTol=2e-9;
    cf.AbsTol=1e-13;
    cf.output_dt=5;
    cf.max_step=30;
    dense=linspace(1e-4,.8,5001);
    for si=1:3
        z=load(fullfile(od,sprintf('seed_%d.mat',seeds(si))));
        res=z.res;
        allruns=[allruns,res.runs];
        
        m=res.test_metrics;
        errors{si}=m.trajectory_rmse;
        testrow=struct('seed',seeds(si),'n_test',96,'n_curve_test',numel(m.test_indices),'mean_curve_RMSE',mean(m.trajectory_rmse),'median_curve_RMSE',median(m.trajectory_rmse),'p95_curve_RMSE',prctile(m.trajectory_rmse,95),'worst_point_error',max(m.max_error),'median_max_error',median(m.max_error),'pca_rank',m.pca_rank,'explained_variance',m.explained,'pca_train_RMSE',mean(m.train_reconstruction_rmse),'Xmol_MAE',m.Xmol_mae,'Xmass_MAE',m.Xmass_mae,'actual_feasible',sum(m.actual_feasible),'predicted_feasible',sum(m.predicted_feasible),'false_feasible',m.false_feasible,'false_infeasible',m.false_infeasible);
        if si==1
            testrows=testrow;
        else
            testrows(si)=testrow;
        end
        seedrows=[seedrows,struct('seed',seeds(si),'training_calls',160,'training_feasible',sum(res.training.feasible),'training_ODE_seconds',sum(res.training.seconds),'shared_GP_fit_seconds',res.shared_fit_seconds,'test_prediction_seconds',res.test_prediction_seconds,'seed_wall_seconds',res.seconds)];
        for j=1:numel(res.runs)
            r=res.runs(j);
            [o,secs]=cached_replay(r.best_x);
            replayseconds=replayseconds+secs;
            f=profile(o,grid);
            fD=profile(o,dense);
            targetD=interp1(grid,r.target.F,dense,'pchip');
            rmse=sqrt(trapz(dense,(fD-targetD).^2)/(dense(end)-dense(1)));
            peak=max(abs(fD-targetD));
            feas=o.Xmol(end)>=.8-1e-8 && o.Xmass(end)>=.8-1e-8 && max(o.V)<=c.max_volume+1e-8;
            if ~feas
                rmse=Inf;
            end
            counter=counter+1;
            fine=[fine,struct('seed',seeds(si),'target',r.target.id,'kind',r.target.kind,'method',r.method,'baseline_RMSE',r.baseline.RMSE,'normal_RMSE',r.RMSE,'fine_RMSE',rmse,'fine_max_error',peak,'feasible',feas,'target_met',feas && rmse<=.02 && peak<=.05,'Xmol',o.Xmol(end),'Xmass',o.Xmass(end),'Vmax',max(o.V),'balance',max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]),'best10',r.best10,'best20',r.best20,'best30',r.best30,'fit_seconds',r.fit_seconds,'prediction_seconds',r.prediction_seconds,'extra_ODE_seconds',r.extra_ode_seconds,'online_seconds',r.total_seconds,'extra_calls',30)];
            details=[details,struct('seed',seeds(si),'target',r.target,'method',r.method,'x',r.best_x,'F',f,'denseF',fD,'t',o.t,'Xmol',o.Xmol,'Ftime',o.Finstant)];
            
        end
        % 数据库检索结果单独回放，不计作新增搜索。
        for ti=1:numel(targets)
            r=res.runs((ti-1)*4+1);
            [o,secs]=cached_replay(r.baseline.x);
            replayseconds=replayseconds+secs;
            fD=profile(o,dense);
            targetD=interp1(grid,targets(ti).F,dense,'pchip');
            rmse=sqrt(trapz(dense,(fD-targetD).^2)/(dense(end)-dense(1)));
            peak=max(abs(fD-targetD));
            counter=counter+1;
            fine=[fine,struct('seed',seeds(si),'target',targets(ti).id,'kind',targets(ti).kind,'method','database','baseline_RMSE',r.baseline.RMSE,'normal_RMSE',r.baseline.RMSE,'fine_RMSE',rmse,'fine_max_error',peak,'feasible',true,'target_met',rmse<=.02 && peak<=.05,'Xmol',o.Xmol(end),'Xmass',o.Xmass(end),'Vmax',max(o.V),'balance',max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]),'best10',rmse,'best20',rmse,'best30',rmse,'fit_seconds',0,'prediction_seconds',0,'extra_ODE_seconds',0,'online_seconds',0,'extra_calls',0)];
        end
    end
    R=struct2table(fine);
    writetable(R,fullfile(od,'design_results.csv'));
    writetable(struct2table(testrows),fullfile(od,'holdout_metrics.csv'));
    writetable(struct2table(seedrows),fullfile(od,'seed_costs.csv'));
    methods={'database','random','sqp','scalar_gp','trajectory_gp'};
    kinds={'all','analytic','reachable'};
    summary=struct([]);
    k=0;
    for ki=1:numel(kinds)
        for mi=1:numel(methods)
            sel=strcmp(R.method,methods{mi});
            if ki>1
                sel=sel & strcmp(R.kind,kinds{ki});
            end
            a=R(sel,:);
            k=k+1;
            summary=[summary,struct('group',kinds{ki},'method',methods{mi},'n',height(a),'mean_RMSE',mean(a.fine_RMSE),'median_RMSE',median(a.fine_RMSE),'mean_max_error',mean(a.fine_max_error),'met_count',sum(a.target_met),'feasible_count',sum(a.feasible),'mean_online_seconds',mean(a.online_seconds),'mean_extra_ODE_seconds',mean(a.extra_ODE_seconds),'mean_fit_seconds',mean(a.fit_seconds),'mean_prediction_seconds',mean(a.prediction_seconds))];
        end
    end
    writetable(struct2table(summary),fullfile(od,'method_summary.csv'));
    paired=struct([]);
    k=0;
    for ki=1:numel(kinds)
        for opponent={'database','random','sqp','scalar_gp'}
            a=R(strcmp(R.method,'trajectory_gp'),:);
            b=R(strcmp(R.method,opponent{1}),:);
            a=sortrows(a,{'seed','target'});
            b=sortrows(b,{'seed','target'});
            assert(isequal(a.seed,b.seed) && isequal(a.target,b.target));
            if ki>1
                sel=strcmp(a.kind,kinds{ki});
                a=a(sel,:);
                b=b(sel,:);
            end
            delta=b.fine_RMSE-a.fine_RMSE;
            rel=delta./b.fine_RMSE;
            k=k+1;
            paired=[paired,struct('group',kinds{ki},'opponent',opponent{1},'n',height(a),'wins_1e4',sum(delta>1e-4),'ties_1e4',sum(abs(delta)<=1e-4),'losses_1e4',sum(delta< -1e-4),'mean_RMSE_improvement',mean(delta),'median_relative_improvement',median(rel),'mean_relative_improvement',mean(rel))];
        end
    end
    writetable(struct2table(paired),fullfile(od,'paired_comparisons.csv'));
    H=struct([]);
    for j=1:numel(allruns)
        rr=allruns(j);
        for n=1:190
            dd=rr.D;
            err=NaN;
            if all(isfinite(dd.F(n,:)))
                err=sqrt(((dd.F(n,:)-rr.target.F).^2)*w);
            end
            H=[H,struct('seed',rr.seed,'target',rr.target.id,'method',rr.method,'evaluation',n,'is_shared_training',n<=160,'RMSE',err,'feasible',dd.feasible(n),'ODE_seconds',dd.seconds(n),'x1',dd.X(n,1),'x2',dd.X(n,2),'x3',dd.X(n,3),'x4',dd.X(n,4))];
            
        end
    end
    writetable(struct2table(H),fullfile(od,'all_search_history.csv'));
    load(fullfile(od,'run_elapsed.mat'),'elapsed');
    numerical=struct('max_normal_fine_RMSE_difference',max(abs(R.normal_RMSE-R.fine_RMSE)),'max_balance_error',max(R.balance),'all_final_process_feasible',all(R.feasible),'fine_unique_calls',numel(finecache),'fine_ODE_seconds',replayseconds);
    online_integrations = 0;
    for k = 1:numel(allruns)
        online_integrations = online_integrations + sum(allruns(k).D.pump_ok(161:end));
    end
    cost=struct('training_calls',480,'test_calls',96,'target_generation_calls',attempts,'online_calls',numel(allruns)*30,'online_actual_ODE_calls',online_integrations,'fine_unique_calls',numel(finecache),'total_calls',480+96+attempts+numel(allruns)*30+numel(finecache),'total_actual_ODE_calls',480+96+attempts+online_integrations+numel(finecache),'run_wall_seconds',elapsed,'fine_ODE_seconds',replayseconds,'test_ODE_seconds',sum(test.seconds),'target_ODE_seconds',target_seconds);
    report=struct('date','2026-10-16','test',testrows,'seed_costs',seedrows,'summary',summary,'paired',paired,'numerical',numerical,'cost',cost);
    fid=fopen(fullfile(od,'results_summary.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
    fclose(fid);
    save(fullfile(od,'analysis.mat'),'R','summary','paired','details','grid','dense','targets','report');
    makeplots;
    disp(struct2table(summary));
    disp(struct2table(paired));
    disp(cost);
    function [o,seconds]=cached_replay(x)
        seconds=0;
        for q=1:numel(finecache)
            if max(abs(finecache(q).x-x))<1e-10
                o=finecache(q).o;
                return;
            end
        end
        t=tic;
        [rates,~,ok,~,cc]=feed_decode(x,cf);
        assert(ok);
        o=simulate_reactor(cc,rates);
        seconds=toc(t);
        finecache=[finecache,struct('x',x,'o',o)];
    end
    function makeplots
        set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',11);
        colors=lines(5);
        fig=figure('Visible','off','Color','w','Position',[100,100,1200,820]);
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        nexttile;
        hold on;
        for s=1:3
            plot(sort(errors{s}),'LineWidth',1.5);
        end;
        yline(.02,'k--','Tracking RMSE threshold');
        xlabel('Independent test feed (sorted)');
        ylabel('Trajectory prediction RMSE');
        legend('Seed 101','Seed 202','Seed 303','Location','northwest');
        title('(a) Held-out full trajectories');
        set(gca,'XGrid','on','YGrid','on');
        nexttile;
        hold on;
        for mi=1:5
            a=R(strcmp(R.method,methods{mi}) & strcmp(R.kind,'reachable'),:);
            scatter(mi+linspace(-.12,.12,height(a)),a.fine_RMSE,28,colors(mi,:),'filled');
        end
        xticks(1:5);
        xticklabels({'Database','Random','SQP','Scalar GP','Trajectory GP'});
        xtickangle(20);
        ylabel('Fine-model tracking RMSE');
        title('(b) Four hidden reachable targets, 3 seeds');
        yline(.02,'k--');
        set(gca,'XGrid','on','YGrid','on');
        nexttile;
        a=sortrows(R(strcmp(R.method,'scalar_gp'),:),{'seed','target'});
        b=sortrows(R(strcmp(R.method,'trajectory_gp'),:),{'seed','target'});
        hold on;
        scatter(a.fine_RMSE,b.fine_RMSE,32,'filled');
        lim=max([a.fine_RMSE;b.fine_RMSE])*1.05;
        plot([0 lim],[0 lim],'k--');
        axis equal;
        xlim([0 lim]);
        ylim([0 lim]);
        xlabel('Scalar GP tracking RMSE');
        ylabel('Trajectory GP tracking RMSE');
        title('(c) Paired final designs');
        set(gca,'XGrid','on','YGrid','on');
        nexttile;
        costs=zeros(4,3);
        for mi=2:5
            a=R(strcmp(R.method,methods{mi}),:);
            costs(mi-1,:)=[mean(a.extra_ODE_seconds),mean(a.fit_seconds),mean(a.online_seconds-a.extra_ODE_seconds-a.fit_seconds)];
        end
        bar(costs,'stacked');
        xticklabels({'Random','SQP','Scalar GP','Trajectory GP'});
        ylabel('Mean online seconds per target');
        legend('ODE','GP fit','Other','Location','northwest');
        title('(d) Includes candidate selection and refitting');
        set(gca,'XGrid','on','YGrid','on');
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,'pilot_summary.png'),'Resolution',360);
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,'pilot_summary.pdf'),'ContentType','vector');
        savefig(fig,fullfile(fd,'pilot_summary.fig'));
        close(fig);
        fig=figure('Visible','off','Color','w','Position',[100,100,1200,780]);
        tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        for ti=[3,6,7,8,9,10]
            nexttile;
            hold on;
            plot(grid,targets(ti).F,'k--','LineWidth',1.8);
            for mi=3:5
                targetMatch = false(size(details));
                for k = 1:numel(details)
                    targetMatch(k) = strcmp(details(k).target.id,targets(ti).id);
                end
                ix=find([details.seed]==202 & strcmp({details.method},methods{mi}) & targetMatch);
                plot(grid,details(ix).F,'Color',colors(mi,:),'LineWidth',1.3);
            end
            xlabel('Polymerization progress');
            ylabel('Instantaneous St fraction');
            title(strrep(targets(ti).id,'_',' '));
            set(gca,'XGrid','on','YGrid','on');
            if ti==3
                legend('Target','SQP','Scalar GP','Trajectory GP','Location','best');
            end
        end
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,'pilot_trajectories.png'),'Resolution',360);
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,'pilot_trajectories.pdf'),'ContentType','vector');
        savefig(fig,fullfile(fd,'pilot_trajectories.fig'));
        close(fig);
    end
end
function f=profile(o,grid)
    valid=o.Xmol>1e-8;
    [u,ix]=unique(o.Xmol(valid));
    F=o.Finstant(valid);
    f=interp1(u,F(ix),grid,'linear');
end
