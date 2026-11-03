function run_initial
    root = fileparts(fileparts(mfilename('fullpath')));
    base = fileparts(root);
    restoredefaultpath;
    addpath(fullfile(root,'matlab','core'));
    maxNumCompThreads(1);
    diary(fullfile(root,'outputs','initial_log.txt'));
    for id = ["constant","A16_K8","A16_K2"]
        for seed = [11 22 33]
            for recipe = ["base","mayo"]
                tag = sprintf('%s_%s_s%d',id,recipe,seed);
                if recipe=="base" && id~="constant"
                    old = load(fullfile(base,'raft_initial_composition','historical',['intervention_' tag '.mat']));
                else
                    old = load(fullfile(base,'raft_initial_composition','outputs',[tag '.mat']));
                end
                reference = old.run;
                dest = fullfile(root,'outputs',['initial_' tag '.mat']);
                if isfile(dest)
                    continue
                end
                reused = recipe=="base" && id~="constant";
                if reused
                    gp = load(fullfile(base,'raft_late_feed','outputs',['gp_' tag '.mat']));
                    run = gp.run;
                else
                    c = reference.c;
                    c.shared_initial = reference.X(1:9,:);
                    c.budget = 450;
                    c.sqp_start_budget = 90;
                    fprintf('START %s %s\n',tag,string(datetime('now')));
                    run = search_feed(c,'gp_constrained',seed);
                end
                assert(isequal(run.X(1:9,:),reference.X(1:9,:)));
                assert(max(abs(run.Y(1:9)-reference.Y(1:9)))<1e-10);
                assert(numel(run.Y)==450);
                save(dest,'run','reference','reused','-v7');
                fprintf('DONE %s GP %.9f SQP %.9f seconds %.1f reused %d\n',tag,run.best_score,reference.best_score,run.seconds,reused);
            end
        end
    end
    diary off
end
