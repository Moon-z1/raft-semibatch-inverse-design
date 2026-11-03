function [score,r,o]=evaluate_feed(x,c)
o=[];[rates,q,ok,~,c]=feed_decode(x,c);
r=struct('score',NaN,'feasible',false,'RMSE',NaN,'max_error',NaN,'fraction_within_2pp',NaN,'target_met',false,...
 'Xmol',NaN,'Xmass',NaN,'Vmax',NaN,'peak_pump_mL_s',max(q,[],'all')*1000,'message','');
if ~ok
 score=2+sum(max(q./c.max_pump-1,0),'all')+sum(max(sum(q,1)/c.max_total_flow-1,0));
 r.message='Pump cap exceeded';r.score=score;return
end
try
 o=simulate_reactor(c,rates);r.Xmol=o.Xmol(end);r.Xmass=o.Xmass(end);r.Vmax=max(o.V);
 violation=max(0,c.min_Xmol-r.Xmol)+max(0,c.min_Xmass-r.Xmass)+max(0,r.Vmax/c.max_volume-1);
 if violation>1e-8,score=1+10*violation;r.message='Production or volume constraint';
 else
  valid=o.Xmol>1e-8;[u,ix]=unique(o.Xmol(valid));F=o.Finstant(valid);F=F(ix);
  grid=linspace(1e-4,c.score_Xmol,801);assert(u(end)>=grid(end),'Insufficient production');
  f=interp1(u,F,grid,'linear');assert(all(isfinite(f)),'Composition domain error');
  residual=f-target_curve(grid,c);
  r.RMSE=sqrt(trapz(grid,residual.^2)/(grid(end)-grid(1)));
  r.max_error=max(abs(residual));r.fraction_within_2pp=trapz(grid,double(abs(residual)<=.02))/(grid(end)-grid(1));
  r.target_met=r.RMSE<=c.rmse_tolerance && r.max_error<=c.peak_tolerance;
  score=r.RMSE;r.feasible=true;
 end
catch ex,score=10;r.message=ex.message;
end
r.score=score;
end
