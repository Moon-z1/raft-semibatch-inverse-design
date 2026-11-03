function make_final_figures
% 读取误差预测与候选排序结果，绘制报告中的对照图。
    % 直接从MATLAB归档结果生成报告图表。
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    fd=fullfile(root,'figures');
    set(groot,'defaultAxesFontName','Microsoft YaHei','defaultTextFontName','Microsoft YaHei','defaultAxesFontSize',10);
    P=readtable(fullfile(root,'outputs','same_data_summary.csv'));
    S=readtable(fullfile(root,'outputs','selection_summary.csv'));
    ids={'A10_K8','A16_K8'};
    vs={'mixed_penalty','separated_response'};
    f=figure('Visible','off','Color','w','Position',[100 100 1120 670]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for ti=1:2
        nexttile;
        vals=zeros(3,2);
        for vi=1:2
            for zi=1:3
                zones={'global','neighborhood','elite_quartile'};
                ix=strcmp(P.target,ids{ti}) & strcmp(P.variant,vs{vi}) & strcmp(P.zone,zones{zi});
                vals(zi,vi)=P.Spearman(ix);
            end
        end
        bar(vals);
        ylim([0 1]);
        xticklabels({'全域','最优邻域','优质四分位'});
        ylabel('Spearman排序相关系数');
        title(sprintf('(%c) %s：同数据预测','a'+ti-1,strrep(ids{ti},'_',' / ')));
        grid on;
        if ti==1
            legend('混合罚分','分离响应','Location','northoutside','Orientation','horizontal');
        end
    end
    nexttile;
    hits=zeros(2,2);
    for ti=1:2
        for vi=1:2
            ix=strcmp(S.target,ids{ti}) & strcmp(S.variant,vs{vi});
            hits(ti,vi)=S.mean_top10_hits(ix);
        end
    end
    bar(hits);
    xticklabels({'中幅陡变','大幅陡变'});
    ylabel('推荐10点中真实前10名的个数');
    ylim([0 10]);
    title('(c) 同一候选池的前10名命中');
    grid on;
    nexttile;
    T=readtable(fullfile(root,'outputs','paired_effects.csv'));
    hold on;
    for ti=1:2
        methods={'random','multistart_sqp','gp_ei','gp_mean'};
        for mi=1:4
            ix=strcmp(T.target,ids{ti}) & strcmp(T.opponent,methods{mi});
            x=mi+(ti-1.5)*.18;
            errorbar(x,1e3*T.mean_gain(ix),1e3*(T.mean_gain(ix)-T.CI95_low(ix)),1e3*(T.CI95_high(ix)-T.mean_gain(ix)),'o','LineWidth',1.1,'Color',[.1 .4 .7]+(ti-1)*[.65 -.1 -.55]);
        end
    end
    yline(0,'k--');
    xlim([.5 4.5]);
    xticks(1:4);
    xticklabels({'随机','SQP','原GP-EI','均值贪心'});
    ylabel('对照RMSE − 约束GP RMSE（×10^{-3}）');
    title('(d) 配对均值与95%自助区间');
    grid on;
    finish(f,'final_ai_evidence');
    D=load(fullfile(base,'raft_trajectory_pilot','outputs','analysis.mat'));
    R=D.R;
    f=figure('Visible','off','Color','w','Position',[100 100 1120 720]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    ms={'database','random','sqp','scalar_gp','trajectory_gp'};
    labels={'数据库','随机','SQP','标量GP','轨迹GP'};
    nexttile;
    hold on;
    for mi=1:5
        a=R(strcmp(R.method,ms{mi}) & strcmp(R.kind,'reachable'),:);
        scatter(mi+linspace(-.12,.12,height(a)),a.fine_RMSE,24,'filled');
    end
    xticks(1:5);
    xticklabels(labels);
    ylabel('精细模型轨迹RMSE');
    title('(a) 4条隐藏可达目标 × 3个种子');
    grid on;
    nexttile;
    a=sortrows(R(strcmp(R.method,'scalar_gp'),:),{'seed','target'});
    b=sortrows(R(strcmp(R.method,'trajectory_gp'),:),{'seed','target'});
    hold on;
    for kind={'reachable','analytic'}
        ix=strcmp(a.kind,kind{1});
        scatter(a.fine_RMSE(ix),b.fine_RMSE(ix),28,'filled');
    end
    plot([0 .1],[0 .1],'k--');
    axis equal;
    xlim([0 .1]);
    ylim([0 .1]);
    xlabel('标量GP：最终轨迹RMSE');
    ylabel('轨迹GP：最终轨迹RMSE');
    legend('可达目标','解析目标','等误差线','Location','northwest');
    title('(b) 相同训练数据与追加预算');
    grid on;
    nexttile;
    targetMatch = false(size(D.details));
    for k = 1:numel(D.details)
        targetMatch(k) = strcmp(D.details(k).target.kind,'reachable');
    end
    ix=find(strcmp({D.details.method},'trajectory_gp') & [D.details.seed]==202 & targetMatch);
    d=D.details(ix(1));
    hold on;
    plot(D.grid,d.target.F,'k--','LineWidth',1.8);
    colors=lines(3);
    for mi=3:5
        targetMatch = false(size(D.details));
        for k = 1:numel(D.details)
            targetMatch(k) = strcmp(D.details(k).target.id,d.target.id);
        end
        j=find(strcmp({D.details.method},ms{mi}) & [D.details.seed]==202 & targetMatch);
        plot(D.grid,D.details(j).F,'Color',colors(mi-2,:),'LineWidth',1.2);
    end
    xlabel('聚合进度 ξ');
    ylabel('瞬时St摩尔分数');
    title('(c) 预定示例：首条可达目标，种子202');
    legend('目标','SQP','标量GP','轨迹GP','Location','best');
    grid on;
    nexttile;
    vals=zeros(4,3);
    for mi=2:5
        a=R(strcmp(R.method,ms{mi}),:);
        vals(mi-1,:)=[mean(a.extra_ODE_seconds),mean(a.fit_seconds),mean(a.online_seconds-a.extra_ODE_seconds-a.fit_seconds)];
    end
    bar(vals,'stacked');
    xticklabels(labels(2:5));
    ylabel('平均在线耗时 / s');
    title('(d) 每目标追加30次评价的完整耗时');
    legend('动力学计算','GP拟合','选点等开销','Location','northwest');
    grid on;
    finish(f,'final_trajectory');
    % 导出每次搜索的分段进料数值。
    addpath(fullfile(base,'raft_trajectory_pilot','matlab','core'),'-begin');
    z=load(fullfile(base,'raft_trajectory_pilot','outputs','protocol_data.mat'),'c');
    rows=struct([]);
    for j=1:numel(D.details)
        d=D.details(j);
        [rates,~,ok,~,c]=feed_decode(d.x,z.c);
        assert(ok);
        for k=1:size(rates,2)
            rows=[rows,struct('target',d.target.id,'method',d.method,'seed',d.seed,'segment',k,'start_s',(k-1)*2400,'end_s',k*2400,'St_mol_s',rates(1,k),'GMA_mol_s',rates(2,k),'x1',d.x(1),'x2',d.x(2),'x3',d.x(3),'x4',d.x(4))];
            
        end
    end
    writetable(struct2table(rows),fullfile(root,'outputs','pilot_best_feeds.csv'));
    function finish(fig,name)
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,[name '.png']),'Resolution',360);
        % 导出前统一坐标框、字体、图例和单位。
        packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(packageFolder);
        apply_report_style(fig);
        exportgraphics(fig,fullfile(fd,[name '.pdf']),'ContentType','vector');
        savefig(fig,fullfile(fd,[name '.fig']));
        close(fig);
    end
end
