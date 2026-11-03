function run_final_evidence
% 在相同数据上比较预测误差、候选排序和实际推荐效果。
    root=fileparts(fileparts(mfilename('fullpath')));
    base=fileparts(root);
    old=fullfile(base,'raft_course_report','outputs');
    od=fullfile(root,'outputs');
    addpath(fullfile(base,'raft_course_report','matlab'),'-begin');
    maxNumCompThreads(1);
    wall=tic;
    ids={'A10_K8','A16_K8'};
    seeds=11:11:110;
    methods={'random','multistart_sqp','gp_ei','gp_mean'};
    paired=struct([]);
    prediction=struct([]);
    selection=struct([]);
    confusion=struct([]);
    allpred=struct([]);
    sample_calls=0;
    for ti=1:2
        id=ids{ti};
        runs=cell(1,10);
        constrained=cell(1,10);
        allX=[];
        anchor=[];
        anchorScore=Inf;
        for si=1:10
            z=load(fullfile(old,sprintf('scan_%s_gp_ei_s%d.mat',id,seeds(si))));
            runs{si}=z.run;
            allX=[allX;z.run.X];
            
            if z.run.best_score<anchorScore
                anchor=z.run.best_x;
                anchorScore=z.run.best_score;
            end
            z=load(fullfile(old,sprintf('scan_%s_gp_constrained_s%d.mat',id,seeds(si))));
            constrained{si}=z.run;
        end
        for mi=1:numel(methods)
            delta=zeros(10,1);
            oldSeconds=delta;
            newSeconds=delta;
            for si=1:10
                z=load(fullfile(old,sprintf('scan_%s_%s_s%d.mat',id,methods{mi},seeds(si))));
                assert(max(abs(z.run.X(1:9,:)-constrained{si}.X(1:9,:)),[],'all')<1e-12);
                delta(si)=z.run.best_score-constrained{si}.best_score;
                oldSeconds(si)=z.run.seconds;
                newSeconds(si)=constrained{si}.seconds;
            end
            rng(2026101700+ti*10+mi);
            boot=mean(delta(randi(10,10,10000)),1);
            ci=prctile(boot,[2.5,97.5]);
            paired=[paired,struct('target',id,'opponent',methods{mi},'n_seeds',10,'mean_gain',mean(delta),'CI95_low',ci(1),'CI95_high',ci(2),'median_gain',median(delta),'wins',sum(delta>0),'historical_opponent_seconds',mean(oldSeconds),'historical_constrained_seconds',mean(newSeconds))];
            
        end
        c=runs{1}.c;
        file=fullfile(od,['fresh_test_' id '.mat']);
        if exist(file,'file')
            d=load(file);
            X=d.X;
            Y=d.Y;
            records=d.records;
            seconds=d.seconds;
        else
            rng(2026101740+ti);
            X=sample_controls(64,c);
            while size(X,1)<128
                x=max(-1,min(1,anchor+.12*randn(1,4)));
                [~,~,ok]=feed_decode(x,c);
                if ok && min(pdist2(x,[allX;X]))>1e-6
                    X=[X;x];
                end 
            end
            assert(min(pdist2(X,allX),[],'all')>1e-6);
            Y=zeros(128,1);
            records=struct([]);
            timer=tic;
            for j=1:128
                [Y(j),rr]=evaluate_feed(X(j,:),c);
                records=[records;rr];
            end 
            seconds=toc(timer);
            save(file,'X','Y','records','seconds','anchor','anchorScore','c');
            sample_calls=sample_calls+128;
        end
        truth=[records.feasible].';
        feasidx=find(truth);
        [~,ix]=sort(Y(feasidx));
        top=feasidx(ix(1:10));
        elite=feasidx(ix(1:ceil(numel(feasidx)*.25)));
        for si=1:10
            r=runs{si};
            rec=r.records;
            valid=isfinite([rec.Xmol].') & isfinite(r.Y);
            feas=valid & [rec.feasible].';
            gp1=fitgp(r.X(valid,:),log(r.Y(valid)),1e-3);
            gp2=fitgp(r.X(feas,:),log(r.Y(feas)),1e-3);
            margins=min([[rec.Xmol].'-c.min_Xmol,[rec.Xmass].'-c.min_Xmass,1-[rec.Vmax].'/c.max_volume],[],2);
            cgp=fitgp(r.X(valid,:),margins(valid),1e-4);
            [cm,cs]=predict(cgp,X);
            predfeas=cm>=0;
            TP=sum(predfeas & truth);
            TN=sum(~predfeas & ~truth);
            FP=sum(predfeas & ~truth);
            FN=sum(~predfeas & truth);
            confusion=[confusion,struct('target',id,'seed',seeds(si),'TP',TP,'TN',TN,'FP',FP,'FN',FN,'actual_infeasible',sum(~truth),'accuracy',(TP+TN)/128,'infeasible_miss_rate',FP/max(1,sum(~truth)))];
            
            for variant=1:2
                labels={'mixed_penalty','separated_response'};
                models={gp1,gp2};
                [mu,sd]=predict(models{variant},X);
                p=exp(mu);
                zones={find(truth & (1:128).'<=64),find(truth & (1:128).'>64),elite};
                zlabels={'global','neighborhood','elite_quartile'};
                for zi=1:3
                    ii=zones{zi};
                    err=p(ii)-Y(ii);
                    rho=NaN;
                    if numel(ii)>2
                        rho=corr(p(ii),Y(ii),'Type','Spearman');
                    end
                    prediction=[prediction,struct('target',id,'seed',seeds(si),'variant',labels{variant},'zone',zlabels{zi},'test_count',numel(ii),'MAE',mean(abs(err)),'RMSE',sqrt(mean(err.^2)),'Spearman',rho,'coverage95',mean(abs(mu(ii)-log(Y(ii)))<=1.96*sd(ii)))];
                    
                end
                candidates=find(predfeas);
                [~,ix]=sort(p(candidates));
                chosen=candidates(ix(1:min(10,numel(ix))));
                best=min(Y(chosen));
                selection=[selection,struct('target',id,'seed',seeds(si),'variant',labels{variant},'recommend_count',numel(chosen),'top10_hits',numel(intersect(chosen,top)),'false_feasible_recommended',sum(~truth(chosen)),'best_recommended_score',best,'gain_vs_training_best',r.best_score-best)];
                
                allpred=[allpred,struct('target',id,'seed',seeds(si),'variant',labels{variant},'pred',p,'logsd',sd,'truth',Y,'feasible',truth,'constraint_mean',cm,'constraint_sd',cs,'chosen',chosen)];
                
            end
        end
        fprintf('SAME-DATA COMPLETE %s, feasible %d/128, fresh ODE seconds %.2f\n',id,sum(truth),seconds);
    end
    writetable(struct2table(paired),fullfile(od,'paired_effects.csv'));
    writetable(struct2table(prediction),fullfile(od,'same_data_prediction.csv'));
    writetable(struct2table(selection),fullfile(od,'top10_selection.csv'));
    writetable(struct2table(confusion),fullfile(od,'constraint_confusion.csv'));
    save(fullfile(od,'same_data_predictions.mat'),'allpred');
    % 串行配对计时，交替方法顺序以减小预热影响。
    warm=study_config;
    [~,~,~]=evaluate_feed(zeros(1,4),warm);
    timing=struct([]);
    for ti=1:2
        for si=1:3
            c=study_config;
            c.target_id=ids{ti};
            c.amplitude=[.10,.16];
            c.amplitude=c.amplitude(ti);
            c.steepness=8;
            c.budget=90;
            order={'multistart_sqp','gp_constrained'};
            if mod(ti+si,2)==0
                order=fliplr(order);
            end
            for mi=1:2
                tag=fullfile(od,sprintf('timing_%s_%s_s%d.mat',ids{ti},order{mi},seeds(si)));
                if exist(tag,'file')
                    d=load(tag);
                    run=d.run;
                else
                    run=search_feed(c,order{mi},seeds(si));
                    save(tag,'run');
                end
                timing=[timing,struct('target',ids{ti},'seed',seeds(si),'method',order{mi},'budget',90,'seconds',run.seconds,'RMSE',run.best_score)];
                
                fprintf('TIMING %s %s s%d %.2fs\n',ids{ti},order{mi},seeds(si),run.seconds);
            end
        end
    end
    writetable(struct2table(timing),fullfile(od,'matched_timing.csv'));
    aggregate_final(root);
    elapsed=toc(wall);
    save(fullfile(od,'supplement_elapsed.mat'),'elapsed','sample_calls');
end
function gp=fitgp(X,y,sigma)
    [xx,ii]=unique(X,'rows','stable');
    gp=fitrgp(xx,y(ii),'KernelFunction','ardmatern52','Standardize',true,'Sigma',sigma,'ConstantSigma',true,'SigmaLowerBound',sigma/1000,'FitMethod','exact','PredictMethod','exact');
end
function aggregate_final(root)
    od=fullfile(root,'outputs');
    P=readtable(fullfile(od,'same_data_prediction.csv'));
    S=readtable(fullfile(od,'top10_selection.csv'));
    C=readtable(fullfile(od,'constraint_confusion.csv'));
    T=readtable(fullfile(od,'matched_timing.csv'));
    agg=struct([]);
    sel=struct([]);
    con=struct([]);
    tim=struct([]);
    for id={'A10_K8','A16_K8'}
        for v={'mixed_penalty','separated_response'}
            for zone={'global','neighborhood','elite_quartile'}
                x=P(strcmp(P.target,id{1}) & strcmp(P.variant,v{1}) & strcmp(P.zone,zone{1}),:);
                agg=[agg,struct('target',id{1},'variant',v{1},'zone',zone{1},'n_models',height(x),'unique_test_points',x.test_count(1),'MAE',mean(x.MAE),'Spearman',mean(x.Spearman,'omitnan'),'coverage95',mean(x.coverage95))];
                
            end
            x=S(strcmp(S.target,id{1}) & strcmp(S.variant,v{1}),:);
            sel=[sel,struct('target',id{1},'variant',v{1},'mean_top10_hits',mean(x.top10_hits),'mean_false_feasible',mean(x.false_feasible_recommended),'mean_best_score',mean(x.best_recommended_score),'mean_actual_gain',mean(x.gain_vs_training_best),'improved_runs_1e4',sum(x.gain_vs_training_best>1e-4))];
            
        end
        x=C(strcmp(C.target,id{1}),:);
        con=[con,struct('target',id{1},'unique_candidates',128,'infeasible_candidates',x.actual_infeasible(1),'mean_TP',mean(x.TP),'mean_TN',mean(x.TN),'mean_FP',mean(x.FP),'mean_FN',mean(x.FN),'mean_accuracy',mean(x.accuracy),'mean_infeasible_miss',mean(x.infeasible_miss_rate),'min_miss',min(x.infeasible_miss_rate),'max_miss',max(x.infeasible_miss_rate))];
        
        for method={'multistart_sqp','gp_constrained'}
            x=T(strcmp(T.target,id{1}) & strcmp(T.method,method{1}),:);
            tim=[tim,struct('target',id{1},'method',method{1},'n_seeds',height(x),'seconds_mean',mean(x.seconds),'seconds_sd',std(x.seconds),'RMSE_mean',mean(x.RMSE))];
            
        end
    end
    writetable(struct2table(agg),fullfile(od,'same_data_summary.csv'));
    writetable(struct2table(sel),fullfile(od,'selection_summary.csv'));
    writetable(struct2table(con),fullfile(od,'confusion_summary.csv'));
    writetable(struct2table(tim),fullfile(od,'timing_summary.csv'));
    report=struct('prediction',agg,'selection',sel,'confusion',con,'timing',tim,'paired',table2struct(readtable(fullfile(od,'paired_effects.csv'))));
    fid=fopen(fullfile(od,'final_evidence.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
    fclose(fid);
end
