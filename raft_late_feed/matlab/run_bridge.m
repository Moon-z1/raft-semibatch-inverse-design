function run_bridge(mode)
    % 复现算法与进料配置的小范围交叉对照。
    if nargin<1
        mode='all';
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    core=fullfile(root,'matlab','core');
    oldpath=path;
    cleanup=onCleanup(@()path(oldpath));
    
    restoredefaultpath;
    addpath(core,'-begin');
    addpath(fullfile(root,'matlab'),'-begin');
    clear search_feed evaluate_feed study_config simulate_reactor feed_decode
    assert(strcmp(which('search_feed'),fullfile(core,'search_feed.m')));
    maxNumCompThreads(1);
    od=fullfile(root,'outputs');
    if any(strcmp(mode,{'all','search'}))
        diary(fullfile(od,'matlab_run_log.txt'));
        closelog=onCleanup(@()diary('off'));
        
        fprintf('Started %s\n',char(datetime('now','TimeZone','Asia/Shanghai')));
        ids={'A16_K8','A16_K2'};
        variants={'base','late_tail'};
        for ti=1:2
            for vi=1:2
                for seed=[11 22 33]
                    tag=sprintf('%s_%s_s%d',ids{ti},variants{vi},seed);
                    dest=fullfile(od,['gp_' tag '.mat']);
                    if exist(dest,'file')
                        fprintf('REUSE completed %s\n',tag);
                        continue;
                    end
                    h=load(fullfile(root,'historical',['intervention_' tag '.mat']));
                    c=h.run.c;
                    assert(c.budget==450 && c.sqp_start_budget==90);
                    global RAFT_BRIDGE_ACCOUNTING RAFT_BRIDGE_PARTIAL
                    RAFT_BRIDGE_ACCOUNTING=struct('score_calls',0,'ode_attempts',0);
                    RAFT_BRIDGE_PARTIAL=[];
                    try
                        run=search_feed(c,'gp_constrained',seed);
                        assert(isequal(run.X(1:9,:),h.run.X(1:9,:)),'Initial controls differ');
                        assert(max(abs(run.Y(1:9)-h.run.Y(1:9)))<1e-10,'Initial scores differ');
                        accounting=RAFT_BRIDGE_ACCOUNTING;
                        assert(accounting.score_calls==451);
                        save(dest,'run','accounting');
                        save_run(fullfile(od,['gp_' tag]),run);
                        % save_run只保存run，需另外补存评价计数记录。
                        save(fullfile(od,['account_' tag '.mat']),'accounting');
                        fprintf('COMPLETE %s RMSE %.9f time %.3f s calls %d ODE %d\n',tag,run.best_score,run.seconds,accounting.score_calls,accounting.ode_attempts);
                    catch ex
                        accounting=RAFT_BRIDGE_ACCOUNTING;
                        partial=RAFT_BRIDGE_PARTIAL;
                        exception=getReport(ex,'extended','hyperlinks','off');
                        save(fullfile(od,['FAILED_' tag '.mat']),'exception','accounting','partial');
                        fprintf(2,'FAILED %s\n%s\n',tag,exception);
                        rethrow(ex)
                    end
                end
            end
        end
    end
    if any(strcmp(mode,{'all','analyze'}))
        analyze_bridge;
    end
    if strcmp(mode,'audit')
        verify_bridge;
    end
end
