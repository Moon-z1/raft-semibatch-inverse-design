function run_scan_constrained
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    restoredefaultpath;
    addpath(fullfile(base,'raft_feed_comparison','matlab','core'));
    if isfolder(fullfile(root,'execution_core'))
        addpath(fullfile(root,'execution_core'));
    end
    maxNumCompThreads(1);
    od=fullfile(root,'outputs');
    diary(fullfile(od,'scan_log.txt'));
    files=dir(fullfile(base,'raft_step45','outputs','scan_*_gp_ei_s*.mat'));
    assert(numel(files)==21);
    for i=1:numel(files)
        tag=strrep(files(i).name,'gp_ei','gp_constrained');
        dest=fullfile(od,tag);
        if isfile(dest)
            continue;
        end
        old=load(fullfile(files(i).folder,files(i).name));
        reference=old.run;
        c=reference.c;
        c.shared_initial=reference.X(1:9,:);
        c.budget=90;
        run=search_feed(c,'gp_constrained',reference.seed);
        assert(isequal(run.X(1:9,:),reference.X(1:9,:)));
        assert(max(abs(run.Y(1:9)-reference.Y(1:9)))<1e-10);
        save(dest,'run','reference','-v7');
        fprintf('DONE %s old %.9f new %.9f\n',tag,reference.best_score,run.best_score);
    end
    diary off
end
