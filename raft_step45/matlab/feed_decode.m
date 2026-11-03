% 将归一化控制变量换算为分段剂量和流量；总投料保持不变。
function [rates,q,ok,weights]=feed_decode(x,c)
    K=c.segments;
    D=2*(K-1);
    x=x(:);
    assert(numel(x)==D && all(abs(x)<=1+1e-8));
    z=c.logit_bound*[x(1:K-1).',0;x(K:end).',0];
    weights=exp(z-max(z,[],2));
    weights=weights./sum(weights,2);
    rates=c.dose.*weights/(c.feed_time/K);
    q=rates.*(c.p.MW./(1000*c.p.rho_monomer_g_mL));
    ok=all(q<=c.max_pump+1e-12,'all') && all(sum(q,1)<=c.max_total_flow+1e-12);
end
