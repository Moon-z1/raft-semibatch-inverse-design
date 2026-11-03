function run_constraint_ai
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    c=study_config;
    ids={'A10_K8','A16_K8'};
    amps=[.10,.16];
    seeds=11:11:110;
    for si=1:10
        maxNumCompThreads(1);
        for ti=1:2
            cc=c;
            cc.target_id=ids{ti};
            cc.amplitude=amps(ti);
            cc.steepness=8;
            tag=fullfile(od,sprintf('scan_%s_gp_constrained_s%d',cc.target_id,seeds(si)));
            if ~exist([tag '.mat'],'file')
                run=search_feed(cc,'gp_constrained',seeds(si));
                save_run(tag,run);
                fprintf('%s %.6f\n',cc.target_id,run.best_score);
            end
        end
    end
end
