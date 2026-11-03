function p = stgma_parameters(caseID)
% 动力学参数参考 Polymers, 2022, 14: 1448；各参数取值及说明见下文。
% Literature-based effective-exchange RAFT model, NOT the previous adduct model.
% Tenorio-Lopez et al., Polymers 2022,14,1448, doi:10.3390/polym14071448.
% Tables 5,6,7,8,9,10. All RAFT coefficients here are L/(mol s).
% The article's fitted effective coefficients MUST NOT be assigned to
% first-order adduct fragmentation constants in the old 119-state model.
if nargin==0, caseID='70-30'; end
p.species = {'styrene','glycidyl methacrylate'};
p.raft_agent = '2-cyano-2-propyl dodecyl trithiocarbonate CPDT';
p.initiator = '1,1-azobis(cyclohexane carbonitrile) ACHN';
p.temperature_C = 103;
p.temperature_K = p.temperature_C + 273.15;
T = p.temperature_K;
R = 8.31446261815324; % J/(mol K), physical constant, not a fitted parameter
% Table 5. Arrhenius kd activation energy is in J/mol after conversion.
p.MW = [104.15;142.20]; % g/mol
p.kd = 2.2798e19*exp(-166900/(R*T)); % 1/s
p.f = 0.60;
p.reactivity_ratio = [.3724;.6834];
kp11=4.266e7*exp(-3909.61/T);
kp22=5.100e6*exp(-2754/T);
p.kp=[kp11,kp11/p.reactivity_ratio(1);kp22/p.reactivity_ratio(2),kp22];
p.ki=[kp11;kp22]; % article assumption ki_i=kp_ii
ktc1=2.0e10*exp(-1553.01/T);
ktc2=2.57e8*exp(-292/T); % article borrows butyl-acrylate value for GMA
p.ktc=[ktc1,sqrt(ktc1*ktc2);sqrt(ktc1*ktc2),ktc2];
p.ktd=p.ktc; % article assumption
% Tables 7 and 9; quoted uncertainties are preserved, not interpreted as SDs.
p.kaa=[101113;484000];
p.kfa=p.kaa;
p.kaa_quoted_plusminus=[2690;1200];
p.kat=[101113;484000;708224;271517]; % order 1,2,3,4
p.kft=[101113;484000;271517;708224];
p.kat_quoted_plusminus=[2690;1200;1756;7223];
p.kft_quoted_plusminus=[2690;1200;7223;1756];
% C(i,j): active i + dormant j -> dormant i + active j.
% Use the single cross exchange pair once, not twice from Table 1.
p.exchange=[p.kat(1),p.kat(3);p.kat(4),p.kat(2)];
% Table 6, densities in g/mL evaluated at Celsius, NOT Kelvin.
p.rho_monomer_g_mL=[.9193-.000665*p.temperature_C;1.09428-.001041*p.temperature_C];
p.rho_polymer_g_mL=[1.09926-2.65e-4*p.temperature_C;1.08];
p.Tg_monomer_K=[156.15;348.15]; % as reported; do not silently correct GMA
p.Tg_polymer_K=[373.15;347.0];
p.alpha_monomer=[.00062;.0004];
p.alpha_polymer=[.00045;.00048];
p.free_volume_at_Tg=.025;
% Table 8: measured recipe concentrations in mol/L, use as printed.
p.case_ids={'100-00','90-10','85-15','70-30','60-40','00-100'};
p.recipes=[8.6445,0,.0297,.0074;7.6091,.8455,.0304,.0076;...
    7.1370,1.2595,.0307,.0077;5.7546,2.4662,.0316,.0079;...
    4.8669,3.2446,.0322,.0080;0,7.5949,.0360,.0090];
% Table 10: per-recipe fitted parameters; not a time-varying composition law.
p.beta_cases=[.3836,.3581,.3979,.23045,.2925,.00081];
p.beta_quoted_plusminus=[.0023,.0214,.0001,5e-5,5e-5,7e-7];
i=find(strcmp(caseID,p.case_ids),1);
assert(~isempty(i),'Unlisted recipe; no silent interpolation of fitted beta.');
p.case_id=caseID;
p.M0=p.recipes(i,1:2).';p.CTA0=p.recipes(i,3);p.I0=p.recipes(i,4);
p.beta=p.beta_cases(i);
p.model_scope='constant-volume batch, intrinsic rates, low-conversion implementation test';
p.autoacceleration_mode='off_low_conversion_only';
p.stop_mass_conversion=.10; % analyst-imposed pilot cutoff; not paper threshold
p.max_time_s=6*3600; % computational cap only
p.RelTol=1e-8;p.AbsTol=1e-12;
p.experimental_validation=false;
end
