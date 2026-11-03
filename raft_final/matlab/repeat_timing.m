function repeat_timing
    % 计时会话中断后，对完整计时组重新计时。
    % 首次计时保留于matched_timing.csv，不删除记录。
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    addpath(fullfile(base,'raft_course_report','matlab'),'-begin');
    maxNumCompThreads(1);
    od=fullfile(root,'outputs');
    c=study_config;
    evaluate_feed(zeros(1,4),c);
    rows=struct([]);
    ids={'A10_K8','A16_K8'};
    seeds=[11 22 33];
    amps=[.10 .16];
    for ti=1:2
        for si=1:3
            c=study_config;
            c.target_id=ids{ti};
            c.amplitude=amps(ti);
            c.steepness=8;
            c.budget=90;
            order={'multistart_sqp','gp_constrained'};
            if mod(ti+si,2)==0
                order=fliplr(order);
            end
            for mi=1:2
                run=search_feed(c,order{mi},seeds(si));
                save(fullfile(od,sprintf('repeat_timing_%s_%s_s%d.mat',ids{ti},order{mi},seeds(si))),'run');
                rows=[rows,struct('target',ids{ti},'seed',seeds(si),'method',order{mi},'budget',90,'seconds',run.seconds,'RMSE',run.best_score)];
                
                fprintf('REPEAT %s %s s%d %.2fs\n',ids{ti},order{mi},seeds(si),run.seconds);
            end
        end
    end
    writetable(struct2table(rows),fullfile(od,'repeat_timing.csv'));
end
