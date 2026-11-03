function [dy,V,Vf,dc]=reactor_rhs(t,y,c,feed)
[V,Vf]=reactor_properties(y,c);p=c.p;
switch c.dc_mode
 case 'continuous_published'
  % Explicit Eq.33-36 interpretation, no unreported onset threshold fitted.
  dc=exp(-c.beta*(1/Vf-1/c.Vf0));
 case 'off'
  dc=1;
 otherwise,error('Unsupported diffusion mode');
end
p.ktc=p.ktc*dc;p.ktd=p.ktd*dc;
dy=V*stgma_rhs(t,y/V,p);dy(1:2)=dy(1:2)+feed;
end
