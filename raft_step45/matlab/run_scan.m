function run_scan
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    c=study_config;
    targets=target_list;
    methods={'random','multistart_sqp','gp_ei'};
    for si=1:3
        maxNumCompThreads(1);
        seed=c.seeds(si);
        for ti=1:numel(targets)
            cc=c;
            cc.amplitude=targets(ti).amplitude;
            cc.steepness=targets(ti).steepness;
            cc.target_id=targets(ti).id;
            for mi=1:3
                kind=methods{mi};
                tag=fullfile(od,sprintf('scan_%s_%s_s%d',cc.target_id,kind,seed));
                if exist([tag '.mat'],'file')
                    continue;
                end
                run=search_feed(cc,kind,seed);
                save_run(tag,run);
                fprintf('%s %s seed%d: %.6f met%d\n',cc.target_id,kind,seed,run.best_score,run.best_record.target_met);
            end
        end
    end
    disp('ALL SCANS COMPLETE');
end
