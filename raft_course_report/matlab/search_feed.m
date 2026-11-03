% 按评分预算搜索进料；SQP回调共享评价次数，达到上限立即停止。
function run=search_feed(c,method,seed,prior)
% 各方法共用同一评分入口，记录每次评价的进料参数和结果。
    % 延续诊断时，评价次数包含先前已经使用的90次。
    if nargin<4
        prior=[];
    end
    rng(seed,'twister');
    D=2*(c.segments-1);
    starts=0;
    cap=Inf;
    timer=tic;
    if isempty(prior)
        X=[zeros(1,D);sample_controls(8,c)];
        Y=zeros(9,1);
        records=[];
        for i=1:9
            [Y(i),r]=evaluate_feed(X(i,:),c);
            records=[records;r];
        end 
        count=9;
        training=true(9,1);
    else
        X=prior.X;
        Y=prior.Y;
        records=prior.records;
        count=numel(Y);
        training=prior.training;
    end
    if strcmp(method,'multistart_sqp')
        [~,ix]=sort(Y);
        points=X(ix(1:3),:);
        while count<c.budget
            starts=starts+1;
            if starts<=3
                x0=points(starts,:);
            else
                x0=sample_controls(1,c);
            end
            cap=min(c.budget,count+27);
            opts=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',27,...
            'FiniteDifferenceStepSize',2e-3,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'MaxIterations',50);
            try
                fmincon(@objective,x0,[],[],[],[],-ones(1,D),ones(1,D),@pumpconstraints,opts);
            catch ex
                if ~strcmp(ex.identifier,'study:budget')
                    rethrow(ex);
                end
            end
        end
    else
        while count<c.budget
            if strcmp(method,'random')
                x=sample_controls(1,c);
            else
                use=find(training);
                if strcmp(method,'gp_constrained')
                    use=find(training & [records.feasible].');
                end
                % 同一控制只保留第一次记录，避免重复输入影响GP拟合。
                [xu,ia]=unique(X(use,:),'rows','stable');
                yu=log(max(Y(use(ia)),1e-6));
                gp=fitrgp(xu,yu,'KernelFunction','ardmatern52','Standardize',true,'Sigma',1e-3,'ConstantSigma',true,...
                'SigmaLowerBound',1e-6,'FitMethod','exact','PredictMethod','exact');
                [~,ib]=min(Y);
                % 候选由全域随机点和当前最佳进料附近的点组成。
                globalPoints = 2*rand(1536,D)-1;
                localPoints = X(ib,:)+.15*randn(512,D);
                localPoints = max(-1,min(1,localPoints));
                candidates = [globalPoints;localPoints];
                keep=false(size(candidates,1),1);
                for j=1:size(candidates,1)
                    [~,~,keep(j)]=feed_decode(candidates(j,:),c);
                end
                candidates=candidates(keep,:);
                if isempty(candidates)
                    candidates=sample_controls(128,c);
                end
                [mu,sd]=predict(gp,candidates);
                sd=max(sd,1e-10);
                % mu为预测均值，sd为预测标准差；以下计算期望改进EI。
                im=min(yu)-mu-.001;
                z=im./sd;
                ei=im.*normcdf(z)+sd.*normpdf(z);
                if strcmp(method,'gp_constrained')
                    margins=min([[records.Xmol].'-c.min_Xmol,[records.Xmass].'-c.min_Xmass,1-[records.Vmax].'/c.max_volume],[],2);
                    known=find(training & isfinite(margins));
                    [cx,ci]=unique(X(known,:),'rows','stable');
                    cy=margins(known(ci));
                    cg=fitrgp(cx,cy,'KernelFunction','ardmatern52','Standardize',true,'Sigma',1e-4,'ConstantSigma',true,...
                    'SigmaLowerBound',1e-7,'FitMethod','exact','PredictMethod','exact');
                    [cm,cs]=predict(cg,candidates);
                    ei=ei.*normcdf(cm./max(cs,1e-10));
                end
                if strcmp(method,'gp_mean')
                    [~,best]=min(mu);
                else
                    [~,best]=max(ei);
                end
                x=candidates(best,:);
            end
            objective(x);
        end
    end
    assert(count==c.budget);
    [best,ib]=min(Y);
    [~,bestrecord,bestout]=evaluate_feed(X(ib,:),c);
    run=struct('c',c,'method',method,'seed',seed,'X',X,'Y',Y,'records',records,'training',training,...
    'best_x',X(ib,:),'best_score',best,'best_record',bestrecord,'bestout',bestout,'seconds',toc(timer),'starts',starts);
    % 求解器接口：共享主函数中的预算计数和结果数组。
    function value=objective(x)
        if count>=min(c.budget,cap)
            error('study:budget','Budget reached');
        end
        [value,r]=evaluate_feed(x,c);
        count=count+1;
        X(count,:)=x(:).';
        Y(count,1)=value;
        records(count,1)=r;
        training(count,1)=true;
    end
    function [g,h]=pumpconstraints(x)
        [~,q]=feed_decode(x,c);
        g=[reshape(q./c.max_pump-1,[],1);(sum(q,1)/c.max_total_flow-1).'];
        h=[];
    end
end
