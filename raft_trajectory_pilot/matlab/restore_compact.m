function restore_compact
    % 展开压缩记录，使analyze_pilot能按原结构读取。
    root=fileparts(fileparts(mfilename('fullpath')));
    od=fullfile(root,'outputs');
    for seed=[101,202,303]
        path=fullfile(od,sprintf('seed_%d.mat',seed));
        if exist(path,'file')
            continue;
        end
        z=load(fullfile(od,sprintf('seed_%d_compact.mat',seed)));
        res=z.res;
        for j=1:numel(res.runs)
            fields=fieldnames(res.training);
            for fi=1:numel(fields)
                f=fields{fi};
                res.runs(j).D.(f)=[res.training.(f);res.runs(j).D.(f)];
            end
        end
        res=rmfield(res,'archive_format');
        save(path,'res','-v7.3');
    end
end
