function [v,r,o]=segment_score(x,c,edges)
    o=[];
    K=numel(edges)-1;
    w=decode(x,K);
    r=struct('score',NaN,'feasible',false,'RMSE',NaN,'max_error',NaN,'target_met',false,'Xmol',NaN,'Xmass',NaN,'Vmax',NaN,'peak_pump_mL_s',NaN,'message','','ode_attempted',false);
    if min(w,[],'all') < -1e-10
        v=2+sum(max(-w,0),'all');
        r.message='Negative dose fraction';
        r.score=v;
        return
    end
    w=max(w,0);
    w=w./sum(w,2);
    rates=c.dose.*w./diff(edges);
    q=rates.*(c.p.MW./(1000*c.p.rho_monomer_g_mL));
    r.peak_pump_mL_s=max(q,[],'all')*1000;
    if any(q>c.max_pump+1e-12,'all') || any(sum(q,1)>c.max_total_flow+1e-12)
        v=2+sum(max(q./c.max_pump-1,0),'all')+sum(max(sum(q,1)/c.max_total_flow-1,0));
        r.message='Pump cap exceeded';
        r.score=v;
        return
    end
    try
        cc=c;
        cc.feed_edges=edges;
        cc.feed_time=edges(end);
        r.ode_attempted=true;
        o=resolution_simulate(cc,rates);
        r.Xmol=o.Xmol(end);
        r.Xmass=o.Xmass(end);
        r.Vmax=max(o.V);
        violation=max(0,c.min_Xmol-r.Xmol)+max(0,c.min_Xmass-r.Xmass)+max(0,r.Vmax/c.max_volume-1);
        if violation>1e-8
            v=1+10*violation;
            r.message='Production or volume constraint';
        else
            [r.RMSE,r.max_error]=metrics(o,c,801);
            r.feasible=true;
            r.target_met=r.RMSE<=.02 && r.max_error<=.05;
            v=r.RMSE;
        end
    catch ex
        v=10;
        r.message=getReport(ex,'basic','hyperlinks','off');
    end
    r.score=v;
end
function w=decode(x,K)
    t=reshape(x,K-1,2).';
    w=[1-sum(t,2),t];
end
function [rmse,peak,grid,f,residual]=metrics(o,c,N)
    valid=o.Xmol>1e-8;
    [u,ix]=unique(o.Xmol(valid));
    F=o.Finstant(valid);
    F=F(ix);
    grid=linspace(1e-4,c.score_Xmol,N);
    assert(u(end)>=grid(end)-1e-8,'Insufficient production');
    f=interp1(u,F,grid,'linear','extrap');
    residual=f-target_curve(grid,c);
    rmse=sqrt(trapz(grid,residual.^2)/diff(grid([1 end])));
    peak=max(abs(residual));
end
