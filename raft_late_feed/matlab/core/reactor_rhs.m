% 半间歇物料衡算：反应消耗与进料共同决定状态变化。
function [dy,V,Vf,dc]=reactor_rhs(t,y,c,feed)
    [V,Vf]=reactor_properties(y,c);
    p=c.p;
    switch c.dc_mode
    case 'continuous_published'
        % 按文献式33—36计算，不额外拟合未报告的起始阈值。
        dc=exp(-c.beta*(1/Vf-1/c.Vf0));
    case 'off'
        dc=1;
    otherwise
        error('Unsupported diffusion mode');
    end
    p.ktc=p.ktc*dc;
    p.ktd=p.ktd*dc;
    dy=V*stgma_rhs(t,y/V,p);
    dy(1:2)=dy(1:2)+feed;
end
