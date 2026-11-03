function make_figures
% 从已有计算记录读取曲线与误差，统一绘图。
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    od=fullfile(root,'outputs');
    ids=["constant","A04_K2","A04_K8","A10_K2","A10_K8","A16_K2","A16_K8"];
    methods=["random","sqp","gp_constrained"];
    Y=zeros(7,3);
    E=Y;
    baseline=zeros(7,1);
    for i=1:7
        for j=1:3
            values=zeros(3,1);
            for k=1:3
                if j==3
                    folder=od;
                else
                    folder=fullfile(base,'raft_step45','outputs');
                end
                file=fullfile(folder,sprintf('scan_%s_%s_s%d.mat',ids(i),methods(j),11*k));
                if ~isfile(file) && j==2
                    file=strrep(file,'_sqp_','_multistart_sqp_');
                end
                a=load(file);
                values(k)=a.run.best_record.RMSE;
                if j==3
                    baseline(i)=a.reference.Y(1);
                end
            end
            Y(i,j)=mean(values);
            E(i,j)=std(values);
        end
    end
    f=figure('Visible','off','Color','w','Position',[100 100 1100 600]);
    hold on;
    plot(1:7,100*baseline,'k--','DisplayName','恒速进料','LineWidth',1.2);
    labels={'随机搜索','多起点SQP','约束GP'};
    colors=[.5 .5 .5;0 .35 .65;.85 .325 .098];
    for j=1:3
        errorbar(1:7,100*Y(:,j),100*E(:,j),'-o','Color',colors(j,:),'DisplayName',labels{j},'LineWidth',1.4);
    end
    xlim([.7 7.3]);
    xticks(1:7);
    xticklabels({'恒定','小幅缓变','小幅陡变','中幅缓变','中幅陡变','大幅缓变','大幅陡变'});
    ylabel('组成RMSE / %');
    set(gca,'FontName','SimSun','FontSize',12);
    legend('Location','best');
    grid on;
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(od,'fig03.png'),'Resolution',360);
    savefig(f,fullfile(od,'fig03.fig'));
    close(f);
    disp(100*(baseline-Y(:,3))./baseline);
    T=readtable(fullfile(base,'raft_feed_comparison','outputs','all_runs.csv'),'TextType','string');
    f=figure('Visible','off','Color','w','Position',[100 100 1100 800]);
    tiledlayout(2,2,'TileSpacing','compact');
    targets=["A16_K8","A16_K2"];
    vars=["four","ten"];
    labels={'四段','十段'};
    col=[.5 .5 .5;0 .35 .65];
    for i=1:2
        nexttile(i);
        hold on;
        for j=1:2
            z=T(T.block=="segments" & T.target==targets(i) & T.configuration==vars(j) & T.method=="SQP",:);
            [~,ib]=min(z.RMSE);
            tag=sprintf('segments_%s_%s_s%d_SQP',targets(i),vars(j),z.seed(ib));
            folder=fullfile(base,'raft_feed_comparison','outputs');
            a=load(fullfile(folder,[tag '.mat']));
            b=load(fullfile(folder,[tag '_fine.mat']));
            r=a.run;
            if j==1
                u=(b.grid-1e-4)/(.8-1e-4);
                target=.7+.16*tanh(r.c.steepness*(u-.5))/tanh(r.c.steepness/2);
                plot(b.grid,target,'k--','DisplayName','目标','LineWidth',1.5);
            end
            plot(b.grid,b.curve,'Color',col(j,:),'LineWidth',1.5,'DisplayName',labels{j});
        end
        ylim([.4 .95]);
        xlabel('聚合进度');
        ylabel('瞬时St摩尔分数');
        title(sprintf('A = 0.16, κ = %d',r.c.steepness));
        legend('Location','best');
        grid on;
        set(gca,'FontName','SimSun','FontSize',11);
        nexttile(i+2);
        hold on;
        K=numel(r.edges)-1;
        t=reshape(r.best_x,K-1,2).';
        w=[1-sum(t,2),t];
        q=r.c.dose.*w./diff(r.edges).*(r.c.p.MW./r.c.p.rho_monomer_g_mL);
        stairs(r.edges/3600,[q(1,:),q(1,end)],'LineWidth',1.5,'DisplayName','St');
        stairs(r.edges/3600,[q(2,:),q(2,end)],'LineWidth',1.5,'DisplayName','GMA');
        xlabel('时间 / h');
        ylabel('进料流量 / (mL·s^{-1})');
        legend('Location','best');
        grid on;
        set(gca,'FontName','SimSun','FontSize',11);
    end
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(od,'fig06.png'),'Resolution',360);
    savefig(f,fullfile(od,'fig06.fig'));
    close(f);
end
