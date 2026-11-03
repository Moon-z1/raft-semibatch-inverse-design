function run_reachability
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    seeds=[11 22 33];
    for si=1:3
        maxNumCompThreads(1);
        seed=seeds(si);
        for ti=1:2
            data=load(fullfile(od,sprintf('known%d_definition.mat',ti)));
            c=data.c;
            methods={'gp_constrained','multistart_sqp'};
            for mi=1:2
                c.budget=270;
                c.sqp_start_budget=90;
                tag=fullfile(od,sprintf('recovery_known%d_%s_s%d',ti,methods{mi},seed));
                if ~exist([tag '.mat'],'file')
                    run=search_feed(c,methods{mi},seed);
                    assert(~any(max(abs(run.X(1:9,:)-data.xtruth),[],2)<1e-9));
                    save_run(tag,run);
                    fprintf('RECOVERY %s %.7f\n',tag,run.best_score);
                end
            end
        end
    end
    disp('RECOVERY RUNS COMPLETE');
end
