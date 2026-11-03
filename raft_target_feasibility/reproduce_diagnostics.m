function reproduce_diagnostics(mode)
% 核验正文4.3的六条解析目标，不改变归档搜索记录。
if nargin==0
    mode='check';
end
root=fileparts(mfilename('fullpath'));
od=fullfile(root,'outputs');
addpath(root,fullfile(root,'core'));
if strcmp(mode,'replay') || ~isfile(fullfile(od,'expanded_trajectories.mat'))
    % 对十段进料重新作精细积分；不重新搜索。
    audit_expanded;
end
files=dir(fullfile(od,'original_t*_seed*.mat'));
assert(numel(files)==18);
rows={};
for i=1:numel(files)
    s=load(fullfile(od,files(i).name));
    r=s.result;
    assert(r.count==900 && size(r.controls,1)==900);
    for k=1:2
        m=r.fine(k);
        e=m.F-m.target;
        rmse=sqrt(trapz(m.grid,e.^2)/(m.grid(end)-m.grid(1)));
        peak=max(abs(e));
        assert(abs(rmse-m.rmse)<1e-12 && abs(peak-m.peak)<1e-12);
        assert(m.feasible && ~m.passed);
        row=struct('target',r.target.id,'seed',r.seed,'candidate',k,...
            'RMSE_percent',100*rmse,'max_error_percent',100*peak,'passed',m.passed);
        rows{end+1}=row;
    end
end
original=struct2table(vertcat(rows{:}));
s=load(fullfile(od,'expanded_trajectories.mat'));
expanded=s.T;
assert(height(expanded)==6 && all(expanded.process_feasible));
assert(sum(expanded.target_met)==4);
expected=["A07_K3","A07_K6","D07_K3","D13_K6"];
assert(isequal(sort(string(expanded.target(expanded.target_met))),sort(expected(:))));
disp(original(original.candidate==1,:));
disp(expanded(:,{'target','RMSE_percent','max_error_percent','target_met'}));
if strcmp(mode,'tables')
    % 导出两组目标的误差汇总。
    writetable(original,fullfile(od,'original_summary.csv'));
    writetable(expanded,fullfile(od,'expanded_summary.csv'));
end
disp('六条解析目标已核验：原控制范围均未达标，扩展条件四条达标、两条未达标。');
end
