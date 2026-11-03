function reproduce_report(mode)
% 默认核验已有结果；不启动优化、不改变报告中的计时。
% 各节计算入口见 README。
if nargin==0
    mode='check';
end
root=fileparts(mfilename('fullpath'));
oldpath=path;
cleanup=onCleanup(@()path(oldpath));
restoredefaultpath;
addpath(root);
if strcmp(mode,'check')
    % 450次对照包含初配方与分段进料两组，共60条最终记录。
    T=readtable(fullfile(root,'raft_feed_comparison','outputs','all_runs.csv'));
    assert(height(T)==60 && all(T.calls==450));
    assert(sum(T.pass)==0);
    % 扩展验证按8条独立目标汇总，不把24次重复当作24条目标。
    G=readtable(fullfile(root,'raft_next','outputs','per_target_summary.csv'));
    assert(numel(unique(G.target))==8);
    % 七条目标的约束GP扫描每条3次、每次90次。
    folder=fullfile(root,'raft_target_scan','outputs');
    files=dir(fullfile(folder,'scan_*_gp_constrained_s*.mat'));
    assert(numel(files)==21);
    for i=1:numel(files)
        z=load(fullfile(folder,files(i).name));
        assert(numel(z.run.Y)==90);
        assert(isequal(z.run.X(1:9,:),z.reference.X(1:9,:)));
        assert(~z.run.best_record.target_met);
    end
    % 两条模型生成目标，各方法3次，共12次恢复。
    folder=fullfile(root,'raft_revision','outputs');
    files=dir(fullfile(folder,'recovery_*.mat'));
    assert(numel(files)==12);
    for i=1:numel(files)
        z=load(fullfile(folder,files(i).name));
        assert(numel(z.run.Y)==270);
        assert(z.run.best_record.target_met);
    end
    disp('已核对450次记录、八条独立目标、七目标扫描和12次恢复。');
elseif strcmp(mode,'diagnostics')
    addpath(fullfile(root,'raft_target_feasibility'));
    reproduce_diagnostics;
elseif strcmp(mode,'figures')
    % 重新导出报告使用的可编辑MATLAB原图，不重新求解。
    files=dir(fullfile(root,'report_figures','*.fig'));
    for i=1:numel(files)
        f=openfig(fullfile(files(i).folder,files(i).name),'invisible');
        [~,name]=fileparts(files(i).name);
        apply_report_style(f);
        exportgraphics(f,fullfile(root,'report_figures',[name '.png']),'Resolution',360);
        close(f);
    end
    disp('已从归档FIG重新导出图件。');
else
    error('mode使用check、figures或diagnostics。');
end
end
