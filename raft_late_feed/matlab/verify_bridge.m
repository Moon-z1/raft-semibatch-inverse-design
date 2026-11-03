function verify_bridge
    % 只读核验；不运行优化、GP拟合或ODE积分。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    R=readtable(fullfile(od,'matched_results.csv'),'TextType','string');
    delta=[];
    peakdelta=[];
    count=0;
    ode=0;
    maxbalance=0;
    for id=["A16_K8","A16_K2"]
        for variant=["base","late_tail"]
            for seed=[11 22 33]
                tag=sprintf('%s_%s_s%d',id,variant,seed);
                d=load(fullfile(od,['gp_' tag '.mat']));
                f=load(fullfile(od,['fine_' tag '.mat']));
                h=load(fullfile(root,'historical',['intervention_' tag '.mat']));
                a=load(fullfile(od,['account_' tag '.mat']));
                assert(numel(d.run.Y)==450 && numel(d.run.records)==450);
                assert(isequal(d.run.X(1:9,:),h.run.X(1:9,:)) && max(abs(d.run.Y(1:9)-h.run.Y(1:9)))<1e-10);
                assert(isequal(d.run.c,h.run.c));
                o=f.o;
                valid=o.Xmol>1e-8;
                [xx,ix]=unique(o.Xmol(valid));
                ff=o.Finstant(valid);
                ff=ff(ix);
                grid=linspace(1e-4,.8,801);
                residual=interp1(xx,ff,grid,'linear')-target_curve(grid,d.run.c);
                rmse=sqrt(trapz(grid,residual.^2)/diff(grid([1 end])));
                pk=max(abs(residual));
                row=R(R.target==id & R.configuration==variant & R.seed==seed & R.method=="Constrained_GP",:);
                assert(height(row)==1);
                delta(end+1)=abs(rmse-row.RMSE);
                peakdelta(end+1)=abs(pk-row.max_error);
                
                nt=sum(o.c.n0)+sum(o.fed(end,:));
                bal=max(abs(o.y(:,1:2)+o.y(:,21:22)-o.c.n0.'-o.fed),[],'all')/nt;
                maxbalance=max(maxbalance,bal);
                count=count+numel(d.run.Y);
                ode=ode+a.accounting.ode_attempts+f.fine_account.ode_attempts;
            end
        end
    end
    assert(count==5400 && max(delta)<1e-12 && max(peakdelta)<1e-12 && maxbalance<1e-8);
    audit=struct('runs',12,'search_budget_calls',count,'completed_ODE_calls',ode,'max_RMSE_recalculation_difference',max(delta),'max_peak_recalculation_difference',max(peakdelta),'max_component_balance',maxbalance,'matching_initial_controls',true,'matching_all_process_settings',true,'optimization_rerun',false,'ODE_rerun',false,'first_failed_analysis_extra_calls_not_included',true);
    fid=fopen(fullfile(od,'independent_audit.json'),'w');
    fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
    fclose(fid);
    disp(audit)
end
