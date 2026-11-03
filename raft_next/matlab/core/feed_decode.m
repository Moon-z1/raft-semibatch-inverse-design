% 将归一化控制变量换算为分段剂量和流量；总投料保持不变。
function [rates,q,ok,weights,c]=feed_decode(x,c)
    [x,c]=apply_design(x,c);
    K=c.segments;
    D=2*(K-1);
    x=x(:);
    late=[0;0];
    if isfield(c,'late_tail') && c.late_tail
        late=(x(D+(1:2))+1)/2;
        x=x(1:D);
    end
    assert(numel(x)==D && all(abs(x)<=1+1e-8));
    z=c.logit_bound*[x(1:K-1).',0;x(K:end).',0];
    weights=exp(z-max(z,[],2));
    weights=weights./sum(weights,2);
    rates=c.dose.*(1-late).*weights/(c.feed_time/K);
    if isfield(c,'late_tail') && c.late_tail
        rates=[rates,c.dose.*late/(c.tail_end-c.feed_time)];
        c.feed_edges=[linspace(0,c.feed_time,K+1),c.tail_end];
        c.feed_time=c.tail_end;
        c.late_fractions=late;
    end
    q=rates.*(c.p.MW./(1000*c.p.rho_monomer_g_mL));
    ok=all(q<=c.max_pump+1e-12,'all') && all(sum(q,1)<=c.max_total_flow+1e-12);
end
