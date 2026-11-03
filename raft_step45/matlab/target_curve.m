% 根据本阶段目标参数计算希望获得的组成曲线。
function f=target_curve(progress,c)
    % 对称梯度在相同评分区间上的均值为0.70。
    u=(progress-1e-4)/(c.score_Xmol-1e-4);
    if c.amplitude==0
        f=c.target_center+zeros(size(u));
    else
        f=c.target_center+c.amplitude*tanh(c.steepness*(u-.5))/tanh(c.steepness/2);
    end
end
