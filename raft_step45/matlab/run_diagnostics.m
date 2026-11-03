function run_diagnostics
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    base=study_config;
    targets=target_list;
    means=zeros(numel(targets),1);
    for ti=1:numel(targets)
        v=zeros(3,1);
        for si=1:3
            d=load(fullfile(od,sprintf('scan_%s_gp_ei_s%d.mat',targets(ti).id,base.seeds(si))));
            v(si)=d.run.best_score;
        end
        means(ti)=mean(v);
    end
    [~,ix]=sort(means,'descend');
    selected=targets(ix(1:2));
    save(fullfile(od,'selected_diagnostic_targets.mat'),'selected','means');
    for si=1:3
        maxNumCompThreads(1);
        seed=base.seeds(si);
        for ti=1:2
            id=selected(ti).id;
            d=load(fullfile(od,sprintf('scan_%s_gp_ei_s%d.mat',id,seed)));
            original=d.run;
            cases={'budget','pump'};
            for ci=1:2
                name=cases{ci};
                tag=fullfile(od,sprintf('diag_%s_%s_s%d',id,name,seed));
                if exist([tag '.mat'],'file')
                    continue;
                end
                c=original.c;
                c.budget=base.diagnostic_budget;
                prior=original;
                if strcmp(name,'pump')
                    c.max_pump=c.max_pump*base.pump_relaxation;
                    c.max_total_flow=c.max_total_flow*base.pump_relaxation;
                    prior.training=~startsWith(string({prior.records.message}).','Pump cap');
                end
                run=search_feed(c,'gp_ei',seed+1000,prior);
                run.parent_seed=seed;
                run.intervention=name;
                run.base90_score=original.best_score;
                save_run(tag,run);
                fprintf('DIAG %s %s seed%d: %.6f from %.6f\n',id,name,seed,run.best_score,original.best_score);
            end
        end
    end
    disp('ALL DIAGNOSTICS COMPLETE');
end
