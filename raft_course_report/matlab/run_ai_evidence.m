function run_ai_evidence
    % 补充对照采用两个目标、十组配对随机初始化和仅均值选点。
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    old=fullfile(fileparts(root),'raft_step45','outputs');
    c=study_config;
    ids={'A10_K8','A16_K8'};
    amps=[.10,.16];
    methods={'random','multistart_sqp','gp_ei','gp_mean'};
    seeds=11:11:110;
    for si=1:numel(seeds)
        maxNumCompThreads(1);
        seed=seeds(si);
        for ti=1:2
            cc=c;
            cc.target_id=ids{ti};
            cc.amplitude=amps(ti);
            cc.steepness=8;
            for mi=1:4
                tag=sprintf('scan_%s_%s_s%d',cc.target_id,methods{mi},seed);
                dest=fullfile(od,tag);
                if exist([dest '.mat'],'file')
                    continue;
                end
                if exist(fullfile(old,[tag '.mat']),'file')
                    d=load(fullfile(old,[tag '.mat']));
                    run=d.run;
                    save_run(dest,run);
                else
                    run=search_feed(cc,methods{mi},seed);
                    save_run(dest,run);
                end
                fprintf('%s %.6f\n',tag,run.best_score);
            end
        end
    end
end
