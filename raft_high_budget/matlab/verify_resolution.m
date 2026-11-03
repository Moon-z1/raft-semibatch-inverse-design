function verify_resolution
    % 从已保存轨迹独立复算，不新增优化。
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'matlab','core'),'-begin');
    od=fullfile(root,'outputs');
    R=readtable(fullfile(od,'results.csv'),'TextType','string');
    D={};
    assert(height(R)==12 && sum(R.calls)==18000);
    for j=1:height(R)
        tag=sprintf('%s_%s_s%d',R.target(j),R.configuration(j),R.seed(j));
        a=load(fullfile(od,[tag '.mat']));
        q=load(fullfile(od,[tag '_fine.mat']));
        r=a.run;
        o=q.o;
        c=r.c;
        K=numel(r.edges)-1;
        x=reshape(r.best_x,K-1,2).';
        w=[1-sum(x,2),x];
        ratio=w(:,1:2)./w(:,3);
        assert(all(ratio>=exp(-4)-1e-6 & ratio<=exp(4)+1e-6,'all'),'Final early-stage bounds violated');
        assert(min(w,[],'all')>=-1e-10);
        rates=c.dose.*w./diff(r.edges);
        volrate=rates.*(c.p.MW./(1000*c.p.rho_monomer_g_mL));
        assert(all(volrate<=c.max_pump+1e-12,'all') && all(sum(volrate,1)<=c.max_total_flow+1e-12));
        assert(max(abs(sum(rates.*diff(r.edges),2)-c.dose))<1e-9);
        valid=o.Xmol>1e-8;
        [xx,ix]=unique(o.Xmol(valid));
        F=o.Finstant(valid);
        F=F(ix);
        tt=o.t(valid);
        tt=tt(ix);
        g=linspace(1e-4,.8,5001);
        u=(g-1e-4)/(.8-1e-4);
        target=.7+.16*tanh(c.steepness*(u-.5))/tanh(c.steepness/2);
        pred=interp1(xx,F,g,'linear');
        e=pred-target;
        assert(all(isfinite(e)));
        rm=sqrt(sum(diff(g).*(e(1:end-1).^2+e(2:end).^2)/2)/(.8-1e-4));
        [pk,ip]=max(abs(e));
        assert(abs(rm-R.fine_RMSE(j))<1e-12 && abs(pk-R.max_error(j))<1e-12);
        assert(numel(r.Y)==1500 && numel(r.records)==1500);
        assert(abs(min(r.Y(1:450))-R.best450(j))<1e-12);
        t80=interp1(xx,tt,.8);
        tpeak=interp1(xx,tt,g(ip));
        n80=interp1(o.t,o.y(:,1:2),t80);
        fed80=interp1(o.t,o.fed,t80);
        b80=interp1(o.t,o.y(:,21:22),t80);
        early=g<=.1;
        middle=g>.1 & g<=.6;
        late=g>.6;
        D{end+1}=struct('target',R.target(j),'configuration',R.configuration(j),'seed',R.seed(j),'search_met_candidates',sum([r.records.target_met]),'RMSE',rm,'max_error',pk,'peak_progress',g(ip),'peak_time_s',tpeak,'peak_signed_error',e(ip),'first_error',e(1),'last_error',e(end),'early_RMSE',sqrt(trapz(g(early),e(early).^2)/range(g(early))),'middle_RMSE',sqrt(trapz(g(middle),e(middle).^2)/range(g(middle))),'late_RMSE',sqrt(trapz(g(late),e(late).^2)/range(g(late))),'time_at_80pct_s',t80,'St_inventory_fraction_at_80pct',n80(1)/sum(n80),'GMA_feed_fraction_at_80pct',fed80(2)/c.dose(2),'cumulative_St_fraction_at_80pct',b80(1)/sum(b80),'max_early_ratio',max(ratio,[],'all'),'min_early_ratio',min(ratio,[],'all'));
        
    end
    D=struct2table(vertcat(D{:}));
    writetable(D,fullfile(od,'residual_diagnostics.csv'));
    E={};
    for target=["A16_K8","A16_K2"]
        a=R(R.target==target & R.configuration=="four",:);
        b=R(R.target==target & R.configuration=="ten",:);
        E{end+1}=struct('target',target,'relative_RMSE_reduction_pct',100*(mean(a.fine_RMSE)-mean(b.fine_RMSE))/mean(a.fine_RMSE),'absolute_RMSE_reduction',mean(a.fine_RMSE)-mean(b.fine_RMSE),'paired_improvements',sum(b.fine_RMSE<a.fine_RMSE),'four_successes',sum(a.target_met),'ten_successes',sum(b.target_met));
        
    end
    writetable(struct2table(vertcat(E{:})),fullfile(od,'configuration_effect.csv'));
    disp(D(:,{'target','configuration','seed','RMSE','max_error','peak_progress','first_error','last_error','search_met_candidates'}));
    v=struct('all_12_saved_trajectories_recomputed',true,'all_final_doses_pumps_early_bounds_verified',true,'search_calls',18000,'new_model_calls_in_this_audit',0,'all_candidate_met_flags',sum(D.search_met_candidates));
    fid=fopen(fullfile(od,'independent_verification.json'),'w');
    fprintf(fid,'%s',jsonencode(v,PrettyPrint=true));
    fclose(fid);
end
