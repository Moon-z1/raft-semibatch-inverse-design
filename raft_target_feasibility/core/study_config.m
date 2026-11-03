function c=study_config
c=stage3_config;
c.logit_bound=4; % expanded identically BEFORE comparisons; pump cap can now bind
c.amplitude=0;c.steepness=0;c.target_id='constant';
c.target_center=.70;c.rmse_tolerance=.02;c.peak_tolerance=.05;
c.diagnostic_budget=270;c.diagnostic_segments=6;c.pump_relaxation=1.5;
c.study_status='conditional computational benchmark; experimentally unvalidated';
end
