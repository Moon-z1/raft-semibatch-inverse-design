function make_resolution_report_figure
    root=fileparts(fileparts(mfilename('fullpath')));
    h=load(fullfile(root,'outputs','analysis.mat'));
    R=h.R;
    S=h.stored;
    f=figure('Visible','off','Color','w','Position',[100 100 1100 800]);
    tiledlayout(2,2,'TileSpacing','compact');
    ids=["A16_K8","A16_K2"];
    col=[.5 .5 .5;0 .35 .65];
    labels={'Four stages','Ten stages'};
    for ti=1:2
        target=ids(ti);
        nexttile(ti);
        hold on;
        for vi=1:2
            vars=["four","ten"];
            variant=vars(vi);
            z=R(R.target==target & R.configuration==variant,:);
            [~,ib]=min(z.fine_RMSE);
            tag=sprintf('%s_%s_s%d',target,variant,z.seed(ib));
            s=S.(matlab.lang.makeValidName(tag));
            if vi==1
                u=(s.grid-1e-4)/(.8-1e-4);
                t=.7+.16*tanh(s.run.c.steepness*(u-.5))/tanh(s.run.c.steepness/2);
                plot(s.grid,t,'k--','DisplayName','Target','LineWidth',1.5);
            end
            plot(s.grid,s.f,'Color',col(vi,:),'LineWidth',1.5,'DisplayName',labels{vi});
        end
        ylim([.45 .95]);
        xlabel('Polymerization progress');
        ylabel('Instantaneous St fraction');
        title(sprintf('A = 0.16, kappa = %d',s.run.c.steepness));
        legend('Location','best');
        grid on;
        set(gca,'FontName','Times New Roman','FontSize',11);
        nexttile(ti+2);
        hold on;
        z=R(R.target==target & R.configuration=="ten",:);
        [~,ib]=min(z.fine_RMSE);
        tag=sprintf('%s_ten_s%d',target,z.seed(ib));
        s=S.(matlab.lang.makeValidName(tag));
        stairs(s.run.edges/3600,[s.q(1,:),s.q(1,end)],'LineWidth',1.5,'DisplayName','St');
        stairs(s.run.edges/3600,[s.q(2,:),s.q(2,end)],'LineWidth',1.5,'DisplayName','GMA');
        xlabel('Time / h');
        ylabel('Feed rate / (mL s^{-1})');
        legend('Location','best');
        grid on;
        set(gca,'FontName','Times New Roman','FontSize',11);
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','resolution_design.png'),'Resolution',360);
    savefig(f,fullfile(root,'figures','resolution_design.fig'));
    close(f);
end
