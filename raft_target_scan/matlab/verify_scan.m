function verify_scan
    % 核对预算、共同初值和精细积分后的精度判定，不重新优化。
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    files=dir(fullfile(od,'scan_*_gp_constrained_s*.mat'));
    assert(numel(files)==21);
    calls=0;
    max_change=0;
    passed=0;
    for i=1:numel(files)
        a=load(fullfile(od,files(i).name));
        r=a.run;
        assert(size(r.X,1)==90 && numel(r.Y)==90);
        assert(isequal(r.X(1:9,:),a.reference.X(1:9,:)));
        assert(max(abs(r.Y(1:9)-a.reference.Y(1:9)))<1e-10);
        f=load(fullfile(od,['fine_' files(i).name]));
        calls=calls+numel(r.Y);
        passed=passed+f.fine.target_met;
        max_change=max(max_change,abs(r.best_record.RMSE-f.fine.RMSE));
        assert(r.best_record.target_met==f.fine.target_met);
    end
    audit=struct('runs',numel(files),'budget_calls',calls,'best_replays_outside_budget',21,'fine_replays_outside_budget',21,'fine_passed',passed,'max_RMSE_change',max_change,'shared_initials_checked',true);
    assert(calls==1890 && passed==0);
    fid=fopen(fullfile(od,'verification.json'),'w');
    fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
    fclose(fid);
    disp(audit);
end
