function prepare_reachability
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    truth=[.3 -.3 -.25 .25;-.3 .3 .25 -.25];
    rows=[];
    for i=1:2
        c=study_config;
        c.target_id=sprintf('known%d',i);
        fine=c;
        fine.RelTol=2e-9;
        fine.AbsTol=1e-13;
        fine.output_dt=5;
        fine.max_step=30;
        [~,r,o]=evaluate_feed(truth(i,:),fine);
        assert(r.feasible,'Prespecified generating schedule is not feasible');
        valid=o.Xmol>1e-8;
        [xx,ix]=unique(o.Xmol(valid));
        fv=o.Finstant(valid);
        fv=fv(ix);
        c.target_grid=linspace(1e-4,.8,5001);
        c.target_values=interp1(xx,fv,c.target_grid,'pchip');
        [self,sr]=evaluate_feed(truth(i,:),c);
        [baseline,br]=evaluate_feed(zeros(1,4),c);
        rows=[rows;struct('target_id',c.target_id,'known_RMSE',self,'baseline_RMSE',baseline,'baseline_met',br.target_met,...
        'generation_Xmol',r.Xmol,'max_error_known',sr.max_error)];
        
        xtruth=truth(i,:);
        save(fullfile(od,[c.target_id '_definition.mat']),'c','xtruth','o');
    end
    writetable(struct2table(rows),fullfile(od,'known_target_precheck.csv'));
    disp(struct2table(rows));
end
