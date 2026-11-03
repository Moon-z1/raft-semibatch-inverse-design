function save_run(path,run)
    save([path '.mat'],'run');
    t=struct2table(run.records);
    t.evaluation=(1:numel(run.Y)).';
    t.best_so_far=cummin(run.Y);
    t.used_for_GP=run.training;
    for j=1:size(run.X,2)
        t.(['x' num2str(j)])=run.X(:,j);
    end
    writetable(t,[path '_history.csv']);
end
