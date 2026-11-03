function verify_fair
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    T=readtable(fullfile(od,'all_runs.csv'),'TextType','string');
    assert(height(T)==60 && sum(~T.reused)==36);
    assert(all(T.calls==450) && sum(T.calls(~T.reused))==16200);
    maxInitialDifference=0;
    maxDoseError=0;
    exceptionCount=0;
    for target=["A16_K8","A16_K2"]
        for seed=[11 22 33]
            for variant=["four","ten"]
                tag=sprintf('segments_%s_%s_s%d_',target,variant,seed);
                a=load(fullfile(od,[tag 'SQP.mat']));
                b=load(fullfile(od,[tag 'GP.mat']));
                assert(isfield(a,'run') && isfield(b,'run'),'Unfinished reservation');
                assert(isequal(a.run.X(1:3,:),b.run.X(1:3,:)));
                maxInitialDifference=max(maxInitialDifference,max(abs(a.run.Y(1:3)-b.run.Y(1:3))));
                assert(a.run.c.sqp_start_budget==90 && b.run.c.sqp_start_budget==90);
                for j=1:numel(a.run.starts)
                    st=a.run.starts{j};
                    assert(st.last_call-st.first_call+1<=90);
                end
                for result={a.run,b.run}
                    r=result{1};
                    K=numel(r.edges)-1;
                    tail=reshape(r.best_x,K-1,2).';
                    w=[1-sum(tail,2),tail];
                    assert(min(w,[],'all')>=-1e-10);
                    maxDoseError=max(maxDoseError,max(abs(sum(r.c.dose.*w,2)-r.c.dose)));
                end
                exceptionCount=exceptionCount+numel(a.exceptions)+numel(b.exceptions);
            end
            a=load(fullfile(od,sprintf('segments_%s_four_s%d_SQP.mat',target,seed)));
            b=load(fullfile(od,sprintf('segments_%s_ten_s%d_SQP.mat',target,seed)));
            for j=1:3
                tail4=reshape(a.run.initial(j,:),3,2).';
                w4=[1-sum(tail4,2),tail4];
                tail10=reshape(b.run.initial(j,:),9,2).';
                w10=[1-sum(tail10,2),tail10];
                expected=[w4(:,1:3),repmat(w4(:,4)/6,1,6),zeros(2,1)];
                assert(max(abs(w10-expected),[],'all')<1e-12);
            end
        end
    end
    for target=["constant","A16_K8","A16_K2"]
        for seed=[11 22 33]
            total=[];
            for recipe=["base","mayo"]
                z=load(fullfile(od,sprintf('initial_%s_%s_s%d.mat',target,recipe,seed)));
                assert(isequal(z.run.X(1:9,:),z.reference.X(1:9,:)));
                maxInitialDifference=max(maxInitialDifference,max(abs(z.run.Y(1:9)-z.reference.Y(1:9))));
                c=z.run.c;
                if isempty(total)
                    total=c.n0+c.dose;
                else
                    assert(max(abs(total-c.n0-c.dose))<1e-12);
                end
            end
        end
    end
    assert(maxInitialDifference<1e-10 && maxDoseError<1e-10);
    S=struct('complete_results',60,'new_runs',36,'new_budget_calls',16200,'new_search_ODE_attempts',sum(T.search_ODE(~T.reused)),'new_best_replays_outside_budget',12,'fine_replays',60,'historical_prefix_replays',12,'max_initial_score_difference',maxInitialDifference,'max_dose_error_mol',maxDoseError,'max_conservation_residual',max(T.balance),'max_search_to_fine_difference',max(T.fine_change),'process_feasible',sum(T.feasible),'composition_passes',sum(T.pass),'segment_exceptions',exceptionCount,'runtime_comparison','Descriptive only: concurrent independent workers and reused historical records');
    S.total_new_ODE_attempts=S.new_search_ODE_attempts+12+60+12;
    S.ODE_count_definition = 'One attempt means one complete feed-scheme forward simulation. Internal ode15s restarts at feed boundaries are not counted separately.';
    % 中断记录不作为额外重复；其计算消耗单列，检查重跑是否一致。
    interrupted = dir(fullfile(od,'interrupted_records','*_partial.mat'));
    savedCalls = 0;
    savedODE = 0;
    maxRestartDifference = 0;
    for j = 1:numel(interrupted)
        prior = load(fullfile(interrupted(j).folder,interrupted(j).name));
        name = strrep(interrupted(j).name,'_partial.mat','.mat');
        current = load(fullfile(od,name));
        n = numel(prior.partial.Y);
        savedCalls = savedCalls + n;
        savedODE = savedODE + prior.partial.odes;
        assert(isequal(prior.partial.X,current.run.X(1:n,:)));
        maxRestartDifference = max(maxRestartDifference,max(abs(prior.partial.Y-current.run.Y(1:n))));
    end
    assert(maxRestartDifference<1e-10);
    S.interrupted_saved_calls = savedCalls;
    S.interrupted_saved_ODE_attempts = savedODE;
    S.all_new_calls_lower_bound = S.new_budget_calls + savedCalls;
    S.all_new_ODE_lower_bound = S.total_new_ODE_attempts + savedODE;
    S.restart_max_score_difference = maxRestartDifference;
    S.interruption_note = 'Incomplete processes exited without final records. Calls after the last checkpoint are unknown; actual total cost is at least this lower bound.';
    fid=fopen(fullfile(od,'verification.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(S,PrettyPrint=true));
    fclose(fid);
    disp(S);
end
