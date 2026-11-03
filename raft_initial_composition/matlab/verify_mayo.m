% 读取已保存记录进行核验和绘图，不重新运行优化。
function verify_mayo
% 检查配方变化前后的总投料和初始方案是否一致。
    % 读取已有记录核验并绘图，不进行搜索或GP拟合。
    root=fileparts(fileparts(mfilename('fullpath')));
    restoredefaultpath;
    addpath(fullfile(root,'matlab/core'),fullfile(root,'matlab'));
    maxNumCompThreads(1);
    s=load(fullfile(root,'outputs','mayo_analysis.mat'));
    T=s.T;
    A=s.A;
    assert(height(T)==18 && all(T.calls==450));
    new=T.method=="mayo"|T.target=="constant";
    calls=sum(T.calls(new));
    assert(calls==5400);
    attempt=0;
    pumpreject=0;
    exceptions=0;
    maxdose=0;
    maxflow=0;
    initialdiff=0;
    rows={};
    for j=1:height(T)
        id=T.target(j);
        method=T.method(j);
        seed=T.seed(j);
        tag=sprintf('%s_%s_s%d',id,method,seed);
        if method=="base" && id~="constant"
            path=fullfile(root,'historical',sprintf('intervention_%s_base_s%d.mat',id,seed));
        else
            path=fullfile(root,'outputs',[tag '.mat']);
        end
        q=load(path);
        r=q.run;
        c=r.c;
        [~,flow,ok]=feed_decode(r.best_x,c);
        assert(ok);
        maxflow=max(maxflow,max(sum(flow,1))*1000);
        assert(all(c.n0>0)&&all(c.dose>0));
        p=stage3_config;
        maxdose=max(maxdose,max(abs(c.n0+c.dose-p.n0-p.dose)));
        assert(c.segments==3 && c.feed_time==7200 && c.budget==450 && c.sqp_start_budget==90);
        if method=="mayo"
            initialdiff=max(initialdiff,max(abs(r.X(1:9,:)-c.shared_initial),[],'all'));
        end
        msg=string({r.records.message});
        reject=sum(msg=="Pump cap exceeded");
        other=msg~="" & msg~="Pump cap exceeded" & msg~="Production or volume constraint";
        rows{end+1}=struct('target',id,'method',method,'seed',seed,'score_calls',numel(r.Y),'ODE_attempts_in_budget',numel(r.Y)-reject,'pump_rejections',reject,'other_errors',sum(other),'new_run',new(j));
        if new(j)
            attempt=attempt+numel(r.Y)-reject;
            pumpreject=pumpreject+reject;
            exceptions=exceptions+sum(other);
        end
    end
    C=struct2table(vertcat(rows{:}));
    writetable(C,fullfile(root,'outputs','evaluation_accounting.csv'));
    fid=fopen(fullfile(root,'outputs','verification.txt'),'w');
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'New optimization runs: 12\nReused baseline runs: 6\nNew in-budget score calls: %d\nODE attempts within new budgets: %d\nPump-rejected scores: %d\nOther errors in new records: %d\nAdditional best replays: 12\nFine replays: 18\n',calls,attempt,pumpreject,exceptions);
    fprintf(fid,'Max fixed-total discrepancy / mol: %.12g\nMax total pump / mL s-1: %.12g\nMax conservation residual: %.12g\nMax fine RMSE change: %.12g\nMax shared-initial difference: %.12g\nMin endpoint mol conversion: %.12g\nMin endpoint mass conversion: %.12g\nMax volume / L: %.12g\n',maxdose,maxflow,max(T.balance),max(abs(T.fine_change)),initialdiff,min(T.Xmol),min(T.Xmass),max(T.Vmax));
    for id=["constant","A16_K8","A16_K2"]
        b=A(A.target==id&A.method=="base",:);
        m=A(A.target==id&A.method=="mayo",:);
        fprintf(fid,'%s mean RMSE decrease / percent: %.9f\n',id,100*(1-m.mean_RMSE/b.mean_RMSE));
    end
    f=figure('Visible','off','Color','w','Position',[100 100 1200 425]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    ids=["constant","A16_K8","A16_K2"];
    names=["(a) Constant","(b) Steep gradient","(c) Gradual gradient"];
    for j=1:3
        nexttile;
        hold on;
        v=s.saved.(sprintf('%s_base_s11',ids(j)));
        c=study_config;
        if j>1
            c.amplitude=.16;
            c.steepness=8;
            if j==3
                c.steepness=2;
            end;
        end
        plot(v.grid,target_curve(v.grid,c)*100,'k--','LineWidth',1.7,'DisplayName','Target');
        for method=["base","mayo"]
            t=T(T.target==ids(j)&T.method==method,:);
            [~,i]=min(t.RMSE);
            v=s.saved.(sprintf('%s_%s_s%d',ids(j),method,t.seed(i)));
            label='Fixed initial recipe';
            color=[.45 .45 .45];
            if method=="mayo"
                label='Mayo-Lewis initial recipe';
                color=[0 .447 .741];
            end
            plot(v.grid,v.F*100,'Color',color,'LineWidth',1.7,'DisplayName',label);
        end
        title(names(j),'FontWeight','normal');
        xlabel('Polymerization progress');
        ylabel('Instantaneous St fraction / %');
        xlim([0 .8]);
        ylim([40 90]);
        set(gca,'XGrid','on','YGrid','on','FontName','Times New Roman','FontSize',13,'Box','on');
        if j==1
            legend('Location','northwest','FontSize',11);
        end
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','mayo_comparison.png'),'Resolution',360);
    savefig(f,fullfile(root,'figures','mayo_comparison.fig'));
    close(f);
    disp(fileread(fullfile(root,'outputs','verification.txt')));
end
