% 固定总投料，对比原初始配方与Mayo–Lewis反算配方。
function run_mayo
% 先确定两组初始配方，再在相同预算下优化进料。
    root=fileparts(fileparts(mfilename('fullpath')));
    restoredefaultpath;
    addpath(fullfile(root,'matlab/core'),fullfile(root,'matlab'));
    maxNumCompThreads(1);
    rows={};
    saved=struct;
    for id=["constant","A16_K8","A16_K2"]
        for seed=[11 22 33]
            if id=="constant"
                tag=sprintf('%s_base_s%d',id,seed);
                dest=fullfile(root,'outputs',[tag '.mat']);
                if exist(dest,'file')
                    s=load(dest);
                    base=s.run;
                else
                    c=study_config;
                    c.budget=450;
                    c.sqp_start_budget=90;
                    h=load(fullfile(root,'historical',sprintf('scan_constant_gp_ei_s%d.mat',seed)));
                    c.anchor_x=h.run.best_x;
                    base=search_feed(c,'multistart_sqp',seed);
                    run=base;
                    save(dest,'run');
                end
            else
                h=load(fullfile(root,'historical',sprintf('intervention_%s_base_s%d.mat',id,seed)));
                base=h.run;
            end
            c=base.c;
            F0=target_curve(1e-4,c);
            r=c.p.reactivity_ratio;
            terminal=@(f)(r(1)*f^2+f*(1-f))/(r(1)*f^2+2*f*(1-f)+r(2)*(1-f)^2);
            f=fzero(@(f)terminal(f)-F0,[.001 .999]);
            total=c.n0+c.dose;
            vm=c.p.MW./(1000*c.p.rho_monomer_g_mL);
            frac=[f;1-f];
            c.n0=frac*c.V0/(frac.'*vm);
            c.dose=total-c.n0;
            c.shared_initial=base.X(1:9,:);
            c.mayo_fraction=f;
            assert(all(c.dose>0)&&max(abs(c.n0+c.dose-total))<1e-12);
            tag=sprintf('%s_mayo_s%d',id,seed);
            dest=fullfile(root,'outputs',[tag '.mat']);
            if exist(dest,'file')
                s=load(dest);
                mayo=s.run;
            else
                fprintf('START %s %s\n',tag,string(datetime('now')));
                mayo=search_feed(c,'multistart_sqp',seed);
                run=mayo;
                save(dest,'run');
            end
            assert(isequal(base.X(1:9,:),mayo.X(1:9,:)));
            for method=["base","mayo"]
                if method=="base"
                    run=base;
                else
                    run=mayo;
                end
                c=run.c;
                c.RelTol=2e-9;
                c.AbsTol=1e-13;
                c.output_dt=5;
                c.max_step=30;
                [~,rec,o]=evaluate_feed(run.best_x,c);
                assert(rec.feasible);
                valid=o.Xmol>1e-8;
                [x,ix]=unique(o.Xmol(valid));
                ff=o.Finstant(valid);
                ff=ff(ix);
                grid=linspace(.0001,.8,5001);
                F=interp1(x,ff,grid);
                res=F-target_curve(grid,c);
                rm=sqrt(trapz(grid,res.^2)/range(grid));
                pk=max(abs(res));
                bal=max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]);
                [~,q]=feed_decode(run.best_x,c);
                tag=sprintf('%s_%s_s%d',id,method,seed);
                fine=struct('o',o,'grid',grid,'F',F,'residual',res,'RMSE',rm,'max_error',pk,'record',rec);
                save(fullfile(root,'outputs',[tag '_fine.mat']),'fine');
                rows{end+1}=struct('target',id,'method',method,'seed',seed,'initial_St',c.n0(1)/sum(c.n0),'RMSE',rm,'max_error',pk,'Xmol',rec.Xmol,'Xmass',rec.Xmass,'Vmax',rec.Vmax,'calls',numel(run.Y),'seconds',run.seconds,'target_met',rm<=.02&&pk<=.05,'first_F',F(1),'last_F',F(end),'balance',bal,'fine_change',rm-run.best_score,'pump_max',max(q,[],'all')*1000);
                
                saved.(tag)=fine;
                fprintf('RESULT %s %.9f %.9f X %.6f\n',tag,rm,pk,rec.Xmol);
            end
        end
    end
    T=struct2table(vertcat(rows{:}));
    writetable(T,fullfile(root,'outputs','mayo_runs.csv'));
    agg={};
    for id=["constant","A16_K8","A16_K2"]
        for method=["base","mayo"]
            t=T(T.target==id&T.method==method,:);
            agg{end+1}=struct('target',id,'method',method,'mean_RMSE',mean(t.RMSE),'sd_RMSE',std(t.RMSE),'mean_max_error',mean(t.max_error),'mean_Xmol',mean(t.Xmol),'passes',sum(t.target_met),'initial_St',t.initial_St(1),'mean_seconds',mean(t.seconds));
            
        end
    end
    A=struct2table(vertcat(agg{:}));
    writetable(A,fullfile(root,'outputs','mayo_summary.csv'));
    save(fullfile(root,'outputs','mayo_analysis.mat'),'T','A','saved');
    disp(A);
    f=figure('Visible','off','Color','w','Position',[100 100 1100 360]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    ids=["constant","A16_K8","A16_K2"];
    names=["Constant","Steep gradient","Gradual gradient"];
    for j=1:3
        nexttile;
        hold on;
        v=saved.(sprintf('%s_base_s11',ids(j)));
        c=study_config;
        if j>1
            c.amplitude=.16;
            c.steepness=8;
            if j==3
                c.steepness=2;
            end;
        end
        plot(v.grid,target_curve(v.grid,c)*100,'k--','LineWidth',1.6,'DisplayName','Target');
        for method=["base","mayo"]
            t=T(T.target==ids(j)&T.method==method,:);
            [~,i]=min(t.RMSE);
            v=saved.(sprintf('%s_%s_s%d',ids(j),method,t.seed(i)));
            plot(v.grid,v.F*100,'LineWidth',1.5,'DisplayName',method);
        end
        title(names(j));
        xlabel('Polymerization progress');
        ylabel('Instantaneous St fraction / %');
        set(gca,'XGrid','on','YGrid','on','FontName','Times New Roman');
        if j==1
            legend('Location','best');
        end
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','mayo_comparison.png'),'Resolution',360);
    savefig(f,fullfile(root,'figures','mayo_comparison.fig'));
    close(f);
end
