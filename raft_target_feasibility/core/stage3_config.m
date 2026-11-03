function c=stage3_config
% Process specifications below are research design choices, not paper data.
c.p=stgma_parameters('70-30');
c.V0=.25; % L, actual additive monomer volume
c.feed_equivalent_L=.75;
c.total_time=24000;c.feed_time=7200;c.segments=3;
c.max_volume=1.05; % L
c.max_pump=[.00018;.00018];c.max_total_flow=.00025; % L/s
c.min_Xmol=.80;c.min_Xmass=.80;c.score_Xmol=.80;
c.target=.70;c.budget=90;c.seeds=[11 22 33];
c.volume_mode='additive';c.dc_mode='continuous_published';
c.RelTol=2e-7;c.AbsTol=1e-11;
c.output_dt=30;c.max_step=120;
f=[.7;.3];vm=c.p.MW./(1000*c.p.rho_monomer_g_mL);
c.n0=f*c.V0/(f.'*vm);
c.dose=f*c.feed_equivalent_L/(f.'*vm);
% Match the paper's CTA/monomer and initiator/monomer ratios.
c.CTA0=sum(c.n0+c.dose)*c.p.CTA0/sum(c.p.M0);
c.I0=sum(c.n0+c.dose)*c.p.I0/sum(c.p.M0);
c.beta=c.p.beta; % fixed 70-30 value; no instantaneous-composition interpolation
c.status='conditional literature-based semibatch extrapolation, not validated experiment';
end
