function [V,Vf,phi]=reactor_properties(y,c)
% State variables are amounts and amount moments. No solvent or outlet.
p=c.p;vm=p.MW./(1000*p.rho_monomer_g_mL);vp=p.MW./(1000*p.rho_polymer_g_mL);
parts=[max(y(1:2),0).*vm;max(y(21:22),0).*vp];
mixV=sum(parts);assert(mixV>0,'Nonpositive material volume');
phi=parts/mixV;
if strcmp(c.volume_mode,'constant'),V=c.V0;else,V=mixV;end
if sum(y(21:22))>1e-15
 f=y(21:22)/sum(y(21:22));mass=y(21:22).*p.MW;w=mass/sum(mass);
else
 f=c.n0/sum(c.n0);w=c.n0.*p.MW/sum(c.n0.*p.MW);
end
TgP=1/sum(w./p.Tg_polymer_K);
alphaP=sum(f.*p.alpha_polymer); % Table 6 uses molar copolymer fractions F_i
Vf=sum((.025+p.alpha_monomer.*(p.temperature_K-p.Tg_monomer_K)).*phi(1:2))...
  +(.025+alphaP*(p.temperature_K-TgP))*sum(phi(3:4));
assert(Vf>0,'Free-volume closure outside domain');
end
