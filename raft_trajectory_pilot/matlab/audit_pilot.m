function audit_pilot
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    addpath(fullfile(root,'matlab','core'),'-begin');
    p=load(fullfile(od,'protocol_data.mat'));
    h=load(fullfile(od,'hidden_target_truth.mat'));
    A=load(fullfile(od,'analysis.mat'));
    R=A.R;
    seeds=[101,202,303];
    min_test_distance=Inf;
    min_hidden_distance=Inf;
    training_equals=true;
    round_counts=true;
    max_balance=0;
    num_badpump=0;
    firstrows=struct([]);
    k=0;
    allcontrols=[];
    original_methods={'trajectory_gp','scalar_gp','sqp','random'};
    for si=1:3
        z=load(fullfile(od,sprintf('seed_%d.mat',seeds(si))));
        res=z.res;
        min_test_distance=min(min_test_distance,min(pdist2(res.training.X,p.test.X),[],'all'));
        truthX=vertcat(h.truth.x);
        min_hidden_distance=min(min_hidden_distance,min(pdist2(res.training.X,truthX),[],'all'));
        allcontrols=[allcontrols;res.training.X];
        
        for j=1:numel(res.runs)
            r=res.runs(j);
            D=r.D;
            training_equals=training_equals && isequaln(D.X(1:160,:),res.training.X) && isequaln(D.F(1:160,:),res.training.F);
            round_counts=round_counts && size(D.X,1)==190;
            max_balance=max(max_balance,max(D.balance,[],'omitnan'));
            num_badpump=num_badpump+sum(~D.pump_ok(161:end));
            e=D.F-r.target.F;
            rmse=sqrt(e.^2*p.w);
            peak=max(abs(e),[],2);
            success=D.feasible & rmse<=.02 & peak<=.05;
            first=find(success,1);
            if isempty(first)
                first=NaN;
            elseif first<=160
                first=0;
            else
                first=first-160;
            end
            k=k+1;
            firstrows=[firstrows,struct('seed',seeds(si),'target',r.target.id,'kind',r.target.kind,'method',r.method,'first_success_extra_call',first,'best10',r.best10,'best20',r.best20,'best30',r.best30)];
            
        end
    end
    assert(training_equals && round_counts && min_test_distance>1e-7 && min_hidden_distance>1e-7);
    assert(all(R.feasible) && all(R.Xmol>=.8-1e-8) && all(R.Xmass>=.8-1e-8) && all(R.Vmax<=1.05+1e-8) && max_balance<1e-5 && max(abs(R.fine_RMSE-R.normal_RMSE))<5e-5);
    F=struct2table(firstrows);
    writetable(F,fullfile(od,'budget_milestones.csv'));
    milestones=struct([]);
    for kind={'analytic','reachable'}
        for method=original_methods
            r=F(strcmp(F.kind,kind{1}) & strcmp(F.method,method{1}),:);
            milestones=[milestones,struct('kind',kind{1},'method',method{1},'mean_best10',mean(r.best10),'mean_best20',mean(r.best20),'mean_best30',mean(r.best30),'met_before_extra',sum(r.first_success_extra_call==0),'met_within10',sum(r.first_success_extra_call<=10),'met_within20',sum(r.first_success_extra_call<=20),'met_within30',sum(r.first_success_extra_call<=30))];
            
        end
    end
    writetable(struct2table(milestones),fullfile(od,'budget_summary.csv'));
    S=struct2table(A.summary);
    a=S(strcmp(S.group,'reachable') & strcmp(S.method,'trajectory_gp'),:);
    b=S(strcmp(S.group,'reachable') & strcmp(S.method,'scalar_gp'),:);
    aa=S(strcmp(S.group,'analytic') & strcmp(S.method,'trajectory_gp'),:);
    bb=S(strcmp(S.group,'analytic') & strcmp(S.method,'scalar_gp'),:);
    audit=struct('status','passed','unchanged_shared_training',training_equals,'all_online_budgets30',round_counts,'minimum_training_test_distance',min_test_distance,'minimum_training_hidden_control_distance',min_hidden_distance,'maximum_history_balance_error',max_balance,'pump_rejected_online_evaluations',num_badpump,'independent_test_controls',96,'curve_test_controls',sum(all(isfinite(p.test.F),2)),'reachable_mean_error_reduction_fraction',1-a.mean_RMSE/b.mean_RMSE,'analytic_mean_error_reduction_fraction',1-aa.mean_RMSE/bb.mean_RMSE,'reachable_online_time_ratio',a.mean_online_seconds/b.mean_online_seconds,'milestones',milestones);
    fid=fopen(fullfile(od,'integrity_audit.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
    fclose(fid);
    % 压缩记录保留完整轨迹历史，共享训练数据只保存一份。
    for si=1:3
        z=load(fullfile(od,sprintf('seed_%d.mat',seeds(si))));
        res=z.res;
        for j=1:numel(res.runs)
            D=res.runs(j).D;
            fields=fieldnames(D);
            for fi=1:numel(fields)
                D.(fields{fi})=D.(fields{fi})(161:end,:);
            end
            res.runs(j).D=D;
        end
        res.archive_format='training holds shared first 160 rows; each runs(j).D holds only the 30 additional rows';
        save(fullfile(od,sprintf('seed_%d_compact.mat',seeds(si))),'res','-v7.3');
        check=load(fullfile(od,sprintf('seed_%d_compact.mat',seeds(si))));
        assert(isequaln(check.res,res));
        for j=1:numel(res.runs)
            fields=fieldnames(res.training);
            for fi=1:numel(fields)
                f=fields{fi};
                assert(isequaln([res.training.(f);res.runs(j).D.(f)],z.res.runs(j).D.(f)));
            end
        end
    end
    disp(audit);
end
