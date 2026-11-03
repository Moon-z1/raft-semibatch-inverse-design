% 产生归一化进料候选；随机数种子由调用入口设置。
function x=sample_controls(n,c)
    D=2*(c.segments-1);
    x=zeros(n,D);
    i=0;
    while i<n
        z=2*rand(1,D)-1;
        [~,~,ok]=feed_decode(z,c);
        if ok
            i=i+1;
            x(i,:)=z;
        end
    end
end
