function analyze_bridge
% 分别汇总进料时段和搜索方法的影响。
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    old=readtable(fullfile(root,'historical','revision_runs.csv'),'TextType','string');
    rows={};
    checks={};
    feedrows={};
    fine_outputs=struct;
    newcalls=0;
    newodes=0;
    for id=["A16_K8","A16_K2"]
        for variant=["base","late_tail"]
            for seed=[11 22 33]
                tag=sprintf('%s_%s_s%d',id,variant,seed);
                h=load(fullfile(root,'historical',['intervention_' tag '.mat']));
                d=load(fullfile(od,['gp_' tag '.mat']));
                a=load(fullfile(od,['account_' tag '.mat']));
                newcalls=newcalls+a.accounting.score_calls;
                newodes=newodes+a.accounting.ode_attempts;
                prior=old(old.run_id=="intervention_"+tag,:);
                assert(height(prior)==1);
                global RAFT_BRIDGE_ACCOUNTING
                RAFT_BRIDGE_ACCOUNTING=struct('score_calls',0,'ode_attempts',0);
                c=d.run.c;
                fine=c;
                fine.RelTol=2e-9;
                fine.AbsTol=1e-13;
                fine.output_dt=5;
                fine.max_step=30;
                cache=fullfile(od,['fine_' tag '.mat']);
                if exist(cache,'file')
                    cached=load(cache);
                    score=cached.score;
                    r=cached.r;
                    o=cached.o;
                    fine_account=cached.fine_account;
                else
                    [score,r,o]=evaluate_feed(d.run.best_x,fine);
                    fine_account=RAFT_BRIDGE_ACCOUNTING;
                    save(cache,'score','r','o','fine_account');
                end
                assert(r.feasible,'Final feed is infeasible');
                valid=o.Xmol>1e-8;
                [xx,ix]=unique(o.Xmol(valid));
                f=o.Finstant(valid);
                f=f(ix);
                grid=linspace(1e-4,.8,5001);
                er=interp1(xx,f,grid,'linear')-target_curve(grid,c);
                rmse5001=sqrt(trapz(grid,er.^2)/diff(grid([1 end])));
                peak5001=max(abs(er));
                balance=max([o.checks.component o.checks.raft o.checks.units o.checks.chains]);
                assert(balance<1e-5 && abs(score-d.run.best_score)<1e-4);
                newcalls=newcalls+fine_account.score_calls;
                newodes=newodes+fine_account.ode_attempts;
                checks{end+1}=struct('id',tag,'fine_change',abs(score-d.run.best_score),'grid_change',abs(score-rmse5001),'balance',balance,'met_801',r.target_met,'met_5001',rmse5001<=.02 && peak5001<=.05,'initial_controls_equal',isequal(d.run.X(1:9,:),h.run.X(1:9,:)));
                
                fine_outputs.(matlab.lang.makeValidName(tag))=o;
                [~,q,~,~,cc]=feed_decode(d.run.best_x,c);
                if isfield(cc,'feed_edges')
                    edges=cc.feed_edges;
                else
                    edges=linspace(0,cc.feed_time,size(q,2)+1);
                end
                for j=1:size(q,2)
                    feedrows{end+1}=struct('id',tag,'start_s',edges(j),'end_s',edges(j+1),'St_mL_s',1000*q(1,j),'GMA_mL_s',1000*q(2,j));
                end 
                for method=["SQP","Constrained_GP"]
                    if method=="SQP"
                        rr=h.run;
                        rmse=prior.RMSE;
                        pk=prior.max_error;
                        met=prior.target_met;
                        xm=prior.Xmol;
                        xmass=prior.Xmass;
                        vol=prior.Vmax;
                        origin="historical";
                    else
                        rr=d.run;
                        rmse=score;
                        pk=r.max_error;
                        met=r.target_met;
                        xm=r.Xmol;
                        xmass=r.Xmass;
                        vol=r.Vmax;
                        origin="new";
                    end
                    rows{end+1}=struct('target',id,'configuration',variant,'seed',seed,'method',method,'origin',origin,'budget',450,'RMSE',rmse,'max_error',pk,'target_met',met,'Xmol',xm,'Xmass',xmass,'Vmax',vol,'best90',min(rr.Y(1:90)),'best270',min(rr.Y(1:270)),'best450',min(rr.Y),'seconds',rr.seconds);
                    
                end
            end
        end
    end
    R=struct2table(vertcat(rows{:}));
    C=struct2table(vertcat(checks{:}));
    A={};
    P={};
    for id=["A16_K8","A16_K2"]
        for method=["SQP","Constrained_GP"]
            for variant=["base","late_tail"]
                t=R(R.target==id & R.method==method & R.configuration==variant,:);
                A{end+1}=struct('target',id,'method',method,'configuration',variant,'n',height(t),'mean_RMSE',mean(t.RMSE),'std_RMSE',std(t.RMSE),'mean_max_error',mean(t.max_error),'worst_max_error',max(t.max_error),'successes',sum(t.target_met),'mean_seconds',mean(t.seconds));
                
            end
            b=R(R.target==id & R.method==method & R.configuration=="base",:);
            l=R(R.target==id & R.method==method & R.configuration=="late_tail",:);
            P{end+1}=struct('target',id,'method',method,'relative_improvement_pct',100*(mean(b.RMSE)-mean(l.RMSE))/mean(b.RMSE),'paired_improvements',sum(l.RMSE<b.RMSE));
            
            P{end}.mean_absolute_improvement=mean(b.RMSE-l.RMSE);
        end
    end
    A=struct2table(vertcat(A{:}));
    P=struct2table(vertcat(P{:}));
    algorithm_effect={};
    budget_effect={};
    for id=["A16_K8","A16_K2"]
        for variant=["base","late_tail"]
            sq=R(R.target==id & R.configuration==variant & R.method=="SQP",:);
            gp=R(R.target==id & R.configuration==variant & R.method=="Constrained_GP",:);
            algorithm_effect{end+1}=struct('target',id,'configuration',variant,'GP_relative_RMSE_improvement_pct',100*(mean(sq.RMSE)-mean(gp.RMSE))/mean(sq.RMSE),'paired_RMSE_improvements',sum(gp.RMSE<sq.RMSE),'SQP_best_fine_RMSE',min(sq.RMSE),'GP_best_fine_RMSE',min(gp.RMSE),'SQP_mean_max_error',mean(sq.max_error),'GP_mean_max_error',mean(gp.max_error));
            
            for method=["SQP","Constrained_GP"]
                z=R(R.target==id & R.configuration==variant & R.method==method,:);
                budget_effect{end+1}=struct('target',id,'configuration',variant,'method',method,'mean_best90',mean(z.best90),'mean_best270',mean(z.best270),'mean_best450',mean(z.best450),'improvement90_to450_pct',100*(mean(z.best90)-mean(z.best450))/mean(z.best90));
                
            end
        end
    end
    E=struct2table(vertcat(algorithm_effect{:}));
    B=struct2table(vertcat(budget_effect{:}));
    writetable(E,fullfile(od,'algorithm_effect.csv'));
    writetable(B,fullfile(od,'budget_checkpoints.csv'));
    writetable(R,fullfile(od,'matched_results.csv'));
    writetable(A,fullfile(od,'aggregate.csv'));
    writetable(P,fullfile(od,'configuration_effect.csv'));
    writetable(C,fullfile(od,'numerical_checks.csv'));
    writetable(struct2table(vertcat(feedrows{:})),fullfile(od,'new_feed_schedules.csv'));
    S=struct('new_search_runs',12,'budget_calls',5400,'new_scoring_calls_including_replays',newcalls,'new_ode_attempts',newodes,'historical_runs_reused',12,'initial_controls_all_equal',all(C.initial_controls_equal),'max_fine_change',max(C.fine_change),'max_grid_change',max(C.grid_change),'max_balance',max(C.balance),'all_classifications_unchanged',all(C.met_801==C.met_5001));
    S.accounting_scope='Successfully recorded searches and replays; excludes first failed analysis';
    S.failed_analysis_extra_ODE_at_least=1;
    S.full_attempt_total_exact=false;
    fid=fopen(fullfile(od,'verification.json'),'w');
    fprintf(fid,'%s',jsonencode(S,PrettyPrint=true));
    fclose(fid);
    save(fullfile(od,'analysis.mat'),'R','A','P','C','S','E','B','fine_outputs');
    make_bridge_figure(R,root);
    disp(A);
    disp(P);
    disp(S);
end
function make_bridge_figure(R,root)
    f=figure('Visible','off','Color','w','Position',[100 100 1080 420]);
    tiledlayout(1,2,'TileSpacing','compact');
    ids=["A16_K8","A16_K2"];
    titles={'A = 0.16, kappa = 8','A = 0.16, kappa = 2'};
    for ti=1:2
        nexttile;
        hold on
        for method=["SQP","Constrained_GP"]
            for variant=["base","late_tail"]
                q=R(R.target==ids(ti) & R.method==method & R.configuration==variant,:);
                y=mean([q.best90 q.best270 q.best450],1);
                label=strrep(method+" / "+variant,'_',' ');
                plot([90 270 450],y,'-o','LineWidth',1.2,'DisplayName',label);
            end
        end
        xlabel('Scoring calls');
        ylabel('Composition RMSE');
        title(titles{ti});
        grid on;
        legend('Location','best','FontSize',8);
        set(gca,'FontName','Times New Roman');
    end
    savefig(f,fullfile(root,'figures','matched_budget_curves.fig'));
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','matched_budget_curves.png'),'Resolution',360);
    % 导出前统一坐标框、字体、图例和单位。
    packageFolder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(packageFolder);
    apply_report_style(f);
    exportgraphics(f,fullfile(root,'figures','matched_budget_curves.pdf'),'ContentType','vector');
    close(f)
end
