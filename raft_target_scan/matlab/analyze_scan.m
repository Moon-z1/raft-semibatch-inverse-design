function analyze_scan
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    restoredefaultpath;
    addpath(fullfile(base,'raft_feed_comparison','matlab','core'));
    if isfolder(fullfile(root,'execution_core'))
        addpath(fullfile(root,'execution_core'));
    end
    files=dir(fullfile(root,'outputs','scan_*_gp_constrained_s*.mat'));
    assert(numel(files)==21);
    rows={};
    for i=1:numel(files)
        z=load(fullfile(files(i).folder,files(i).name));
        r=z.run;
        c=r.c;
        c.RelTol=2e-9;
        c.AbsTol=1e-13;
        c.output_dt=5;
        c.max_step=30;
        [~,fine,o]=evaluate_feed(r.best_x,c);
        save(fullfile(root,'outputs',['fine_' files(i).name]),'fine','o');
        rows{end+1}=struct('target',string(r.c.target_id),'seed',r.seed,'old_RMSE',z.reference.best_record.RMSE,'RMSE',r.best_record.RMSE,'peak',r.best_record.max_error,'pass',r.best_record.target_met,'fine_RMSE',fine.RMSE,'fine_peak',fine.max_error,'fine_pass',fine.target_met,'seconds',r.seconds);
    end
    T=struct2table(vertcat(rows{:}));
    writetable(T,fullfile(root,'outputs','scan_pairs.csv'));
    rows={};
    ids=unique(T.target,'stable');
    for id=ids.'
        a=T(T.target==id,:);
        rows{end+1}=struct('target',id,'old_mean',mean(a.old_RMSE),'mean_RMSE',mean(a.RMSE),'sd_RMSE',std(a.RMSE),'mean_peak',mean(a.peak),'relative_change_percent',100*(mean(a.RMSE)-mean(a.old_RMSE))/mean(a.old_RMSE),'passes',sum(a.pass),'fine_passes',sum(a.fine_pass),'max_fine_difference',max(abs(a.RMSE-a.fine_RMSE)));
    end
    S=struct2table(vertcat(rows{:}));
    writetable(S,fullfile(root,'outputs','scan_summary.csv'));
    disp(S);
end
