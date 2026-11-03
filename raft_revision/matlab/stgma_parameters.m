% 动力学参数集中于本文件；单位及文献来源见各参数旁的说明。
% 动力学参数参考 Polymers, 2022, 14: 1448；各参数取值及说明见下文。
function p = stgma_parameters(caseID)
    % 采用文献有效交换RAFT模型，与早期加合物模型区分。
    % Tenorio-Lopez et al., Polymers 2022,14,1448, doi:10.3390/polym14071448.
    % 参数见文献表5—10；此处RAFT系数单位均为L/(mol s)。
    % 文献拟合的有效系数不可直接赋给
    % 旧119状态模型的一阶加合物裂解常数。
    if nargin==0
        caseID='70-30';
    end
    p.species = {'styrene','glycidyl methacrylate'};
    p.raft_agent = '2-cyano-2-propyl dodecyl trithiocarbonate CPDT';
    p.initiator = '1,1-azobis(cyclohexane carbonitrile) ACHN';
    p.temperature_C = 103;
    p.temperature_K = p.temperature_C + 273.15;
    T = p.temperature_K;
    R = 8.31446261815324;
    % 单位J/(mol K)，为物理常数，不参与拟合。
    % 文献表5：Arrhenius式中的kd活化能换算为J/mol。
    p.MW = [104.15;142.20];
    % g/mol
    p.kd = 2.2798e19*exp(-166900/(R*T));
    % 1/s
    p.f = 0.60;
    p.reactivity_ratio = [.3724;.6834];
    kp11=4.266e7*exp(-3909.61/T);
    kp22=5.100e6*exp(-2754/T);
    p.kp=[kp11,kp11/p.reactivity_ratio(1);kp22/p.reactivity_ratio(2),kp22];
    p.ki=[kp11;kp22];
    % 沿用文献假设ki_i=kp_ii。
    ktc1=2.0e10*exp(-1553.01/T);
    ktc2=2.57e8*exp(-292/T);
    % 文献将丙烯酸丁酯的数值借用于GMA。
    p.ktc=[ktc1,sqrt(ktc1*ktc2);sqrt(ktc1*ktc2),ktc2];
    p.ktd=p.ktc;
    % 沿用文献假设。
    % 保留文献表7、9的不确定度，不擅自解释为标准差。
    p.kaa=[101113;484000];
    p.kfa=p.kaa;
    p.kaa_quoted_plusminus=[2690;1200];
    p.kat=[101113;484000;708224;271517];
    % 按配方1、2、3、4排列。
    p.kft=[101113;484000;271517;708224];
    p.kat_quoted_plusminus=[2690;1200;1756;7223];
    p.kft_quoted_plusminus=[2690;1200;7223;1756];
    % C(i,j)：活性链i与休眠链j交换，形成休眠链i和活性链j。
    % 交叉交换对只计一次，不能按表1重复计入。
    p.exchange=[p.kat(1),p.kat(3);p.kat(4),p.kat(2)];
    % 文献表6密度单位g/mL；温度代入摄氏度而非开尔文。
    p.rho_monomer_g_mL=[.9193-.000665*p.temperature_C;1.09428-.001041*p.temperature_C];
    p.rho_polymer_g_mL=[1.09926-2.65e-4*p.temperature_C;1.08];
    p.Tg_monomer_K=[156.15;348.15];
    % 保留文献报告值，不自行修改GMA参数。
    p.Tg_polymer_K=[373.15;347.0];
    p.alpha_monomer=[.00062;.0004];
    p.alpha_polymer=[.00045;.00048];
    p.free_volume_at_Tg=.025;
    % 文献表8：配方浓度单位mol/L，沿用原报告数值。
    p.case_ids={'100-00','90-10','85-15','70-30','60-40','00-100'};
    p.recipes=[8.6445,0,.0297,.0074;7.6091,.8455,.0304,.0076;...
    7.1370,1.2595,.0307,.0077;5.7546,2.4662,.0316,.0079;...
    4.8669,3.2446,.0322,.0080;0,7.5949,.0360,.0090];
    % 文献表10为各配方拟合参数，不能解释为随组成变化的定律。
    p.beta_cases=[.3836,.3581,.3979,.23045,.2925,.00081];
    p.beta_quoted_plusminus=[.0023,.0214,.0001,5e-5,5e-5,7e-7];
    i=find(strcmp(caseID,p.case_ids),1);
    assert(~isempty(i),'Unlisted recipe; no silent interpolation of fitted beta.');
    p.case_id=caseID;
    p.M0=p.recipes(i,1:2).';
    p.CTA0=p.recipes(i,3);
    p.I0=p.recipes(i,4);
    p.beta=p.beta_cases(i);
    p.model_scope='constant-volume batch, intrinsic rates, low-conversion implementation test';
    p.autoacceleration_mode='off_low_conversion_only';
    p.stop_mass_conversion=.10;
    % 本研究试验采用的截断值，不是文献阈值。
    p.max_time_s=6*3600;
    % 仅作为计算上限。
    p.RelTol=1e-8;
    p.AbsTol=1e-12;
    p.experimental_validation=false;
end
