function analyze_fair(allow_partial)
% 将各次搜索的最佳进料重新积分，再按目标和配置汇总。
    if nargin<1
        allow_partial=false;
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    restoredefaultpath;
    addpath(fullfile(root,'matlab','core'),fullfile(root,'matlab'));
    maxNumCompThreads(1);
    od = fullfile(root,'outputs');
    rows = {};
    files = dir(fullfile(od,'initial_*.mat'));
    for i = 1:numel(files)
        if contains(files(i).name,'_fine')
            continue
        end
        data = load(fullfile(od,files(i).name));
        tag = erase(files(i).name,'.mat');
        for method = ["SQP","GP"]
            if method=="SQP"
                r = data.reference;
                reused = true;
            else
                r = data.run;
                reused = data.reused;
            end
            tokens = split(tag,'_');
            if tokens(2)=="constant"
                target = "constant";
                recipe = tokens(3);
            else
                target = tokens(2)+"_"+tokens(3);
                recipe = tokens(4);
            end
            finefile = fullfile(od,[tag '_' char(method) '_fine.mat']);
            if isfile(finefile)
                f = load(finefile);
            else
                c = fine_config(r.c);
                [~,rec,o] = evaluate_feed(r.best_x,c);
                f = replay(o,rec,c);
                save(finefile,'-struct','f');
            end
            nODE = sum(~strcmp({r.records.message},'Pump cap exceeded'));
            row = row_data("initial",target,recipe,method,r.seed,r,f,reused,nODE);
            rows{end+1} = row;
        end
    end
    files = dir(fullfile(od,'segments_*.mat'));
    for i = 1:numel(files)
        if contains(files(i).name,{'partial','fine'})
            continue
        end
        data = load(fullfile(od,files(i).name));
        if ~isfield(data,'run') || ~isstruct(data.run)
            assert(isfield(data,'reservation'));
            continue
        end
        r = data.run;
        tag = erase(files(i).name,'.mat');
        finefile = fullfile(od,[tag '_fine.mat']);
        if isfile(finefile)
            f = load(finefile);
        else
            c = fine_config(r.c);
            [~,rec,o] = segment_score(r.best_x,c,r.edges);
            f = replay(o,rec,c);
            save(finefile,'-struct','f');
        end
        row = row_data("segments",r.target,r.variant,r.method,r.seed,r,f,false,r.ode_attempts);
        rows{end+1} = row;
    end
    T = struct2table(vertcat(rows{:}));
    writetable(T,fullfile(od,'all_runs.csv'));
    groups = unique(T(:,{'block','target','configuration','method'}),'rows','stable');
    agg = {};
    for i = 1:height(groups)
        g = groups(i,:);
        z = T(T.block==g.block & T.target==g.target & T.configuration==g.configuration & T.method==g.method,:);
        agg{end+1} = struct('block',g.block,'target',g.target,'configuration',g.configuration,'method',g.method,'runs',height(z),'mean_RMSE',mean(z.RMSE),'sd_RMSE',std(z.RMSE),'mean_peak',mean(z.peak),'passes',sum(z.pass),'feasible',sum(z.feasible),'mean_seconds',mean(z.seconds));
    end
    A = struct2table(vertcat(agg{:}));
    writetable(A,fullfile(od,'summary.csv'));
    save(fullfile(od,'fair_analysis.mat'),'T','A');
    paired={};
    settings=unique(T(:,{'block','target','configuration'}),'rows','stable');
    for i=1:height(settings)
        g=settings(i,:);
        a=T(T.block==g.block & T.target==g.target & T.configuration==g.configuration & T.method=="SQP",:);
        b=T(T.block==g.block & T.target==g.target & T.configuration==g.configuration & T.method=="GP",:);
        [shared,ia,ib]=intersect(a.seed,b.seed);
        if isempty(shared)
            continue
        end
        a=a(ia,:);
        b=b(ib,:);
        delta=a.RMSE-b.RMSE;
        paired{end+1}=struct('block',g.block,'target',g.target,'configuration',g.configuration,'pairs',numel(shared),'SQP_RMSE',mean(a.RMSE),'GP_RMSE',mean(b.RMSE),'GP_relative_gain_percent',100*(mean(a.RMSE)-mean(b.RMSE))/mean(a.RMSE),'GP_lower_runs',sum(delta>0),'SQP_peak',mean(a.peak),'GP_peak',mean(b.peak),'SQP_passes',sum(a.pass),'GP_passes',sum(b.pass));
    end
    if ~isempty(paired)
        P=struct2table(vertcat(paired{:}));
        writetable(P,fullfile(od,'paired_algorithms.csv'));
    end
    effects={};
    for block=["initial","segments"]
        if block=="initial"
            ids=["constant","A16_K8","A16_K2"];
            before="base";
            after="mayo";
        else
            ids=["A16_K8","A16_K2"];
            before="four";
            after="ten";
        end
        for target=ids
            for method=["SQP","GP"]
                a=T(T.block==block & T.target==target & T.method==method & T.configuration==before,:);
                b=T(T.block==block & T.target==target & T.method==method & T.configuration==after,:);
                [shared,ia,ib]=intersect(a.seed,b.seed);
                if isempty(shared)
                    continue
                end
                a=a(ia,:);
                b=b(ib,:);
                effects{end+1}=struct('block',block,'target',target,'method',method,'pairs',numel(shared),'before_RMSE',mean(a.RMSE),'after_RMSE',mean(b.RMSE),'relative_gain_percent',100*(mean(a.RMSE)-mean(b.RMSE))/mean(a.RMSE),'before_peak',mean(a.peak),'after_peak',mean(b.peak));
            end
        end
    end
    if ~isempty(effects)
        writetable(struct2table(vertcat(effects{:})),fullfile(od,'operation_effects.csv'));
    end
    prefix=readtable(fullfile(od,'historical_prefix450.csv'),'TextType','string');
    old=readtable(fullfile(fileparts(root),'raft_high_budget','outputs','aggregate.csv'),'TextType','string');
    budgetRows={};
    for target=["A16_K8","A16_K2"]
        for variant=["four","ten"]
            a=T(T.block=="segments" & T.target==target & T.configuration==variant & T.method=="SQP",:);
            b=prefix(prefix.target==target & prefix.configuration==variant,:);
            c=old(old.target==target & old.configuration==variant,:);
            if height(a)==3
                budgetRows{end+1}=struct('target',target,'configuration',variant,'RMSE_450_start90',mean(a.RMSE),'RMSE_450_start500',mean(b.fine_RMSE),'RMSE_1500_start500',c.mean_RMSE,'peak_450_start90',mean(a.peak),'peak_450_start500',mean(b.peak),'peak_1500_start500',c.mean_max_error);
            end
        end
    end
    if ~isempty(budgetRows)
        writetable(struct2table(vertcat(budgetRows{:})),fullfile(od,'budget_comparison.csv'));
    end
    if ~allow_partial
        assert(height(T)==60,'Missing comparison runs');
        assert(sum(T.calls(~T.reused))==16200);
    end
    assert(all(T.calls==450));
    assert(all(T.balance<1e-8));
    disp(A);
end
function c = fine_config(c)
    c.RelTol=2e-9;
    c.AbsTol=1e-13;
    c.output_dt=5;
    c.max_step=30;
end
function f = replay(o,rec,c)
    assert(~isempty(o),'Fine replay failed');
    valid=o.Xmol>1e-8;
    [u,ix]=unique(o.Xmol(valid));
    F=o.Finstant(valid);
    grid=linspace(.0001,.8,5001);
    curve=interp1(u,F(ix),grid,'linear');
    assert(all(isfinite(curve)));
    error=curve-target_curve(grid,c);
    f=struct('RMSE',sqrt(trapz(grid,error.^2)/range(grid)),'peak',max(abs(error)),'grid',grid,'curve',curve,'rec',rec,'o',o,'balance',max([o.checks.component,o.checks.raft,o.checks.units,o.checks.chains]));
end
function a = row_data(block,target,configuration,method,seed,r,f,reused,nODE)
    a=struct('block',string(block),'target',string(target),'configuration',string(configuration),'method',string(method),'seed',seed,'RMSE',f.RMSE,'peak',f.peak,'feasible',f.rec.feasible,'pass',f.rec.feasible && f.RMSE<=.02 && f.peak<=.05,'Xmol',f.rec.Xmol,'Xmass',f.rec.Xmass,'Vmax',f.rec.Vmax,'seconds',r.seconds,'calls',numel(r.Y),'search_ODE',nODE,'reused',reused,'search_score',r.best_score,'fine_change',abs(f.RMSE-r.best_score),'balance',f.balance);
    if block=="initial"
        nInitial=9;
    else
        nInitial=3;
    end
    a.initial_best=min(r.Y(1:nInitial));
    a.search_improvement=a.initial_best-r.best_score;
    a.process_feasible_calls=sum([r.records.feasible]);
    a.resume_from=0;
    a.resume_seconds=NaN;
    if isfield(r,'resume_from')
        a.resume_from=r.resume_from;
    end
    if isfield(r,'resume_seconds')
        a.resume_seconds=r.resume_seconds;
    end
end
