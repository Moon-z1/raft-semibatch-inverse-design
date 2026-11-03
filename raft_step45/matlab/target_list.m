function targets=target_list
    targets=struct('id','constant','amplitude',0,'steepness',0);
    for a=[.04 .10 .16]
        for k=[2 8]
            targets(end+1)=struct('id',sprintf('A%02d_K%d',round(a*100),k),'amplitude',a,'steepness',k);
            
        end
    end
end
