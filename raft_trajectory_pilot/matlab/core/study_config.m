% 本阶段的目标、进料时间及过程约束；不同阶段的配置不可混用。
function c=study_config
    c=stage3_config;
    c.logit_bound=4;
    % 各方法比较前统一扩展范围，此时泵上限可能成为有效约束。
    c.amplitude=0;
    c.steepness=0;
    c.target_id='constant';
    c.target_center=.70;
    c.rmse_tolerance=.02;
    c.peak_tolerance=.05;
    c.diagnostic_budget=270;
    c.diagnostic_segments=6;
    c.pump_relaxation=1.5;
    c.study_status='conditional computational benchmark; experimentally unvalidated';
end
