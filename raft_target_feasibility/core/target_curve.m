function f=target_curve(progress,c)
if isfield(c,'target_grid'),f=interp1(c.target_grid,c.target_values,progress,'pchip');return;end
% Symmetric ramp has mean .70 on the identical scored production interval.
u=(progress-1e-4)/(c.score_Xmol-1e-4);
if c.amplitude==0,f=c.target_center+zeros(size(u));else
 f=c.target_center+c.amplitude*tanh(c.steepness*(u-.5))/tanh(c.steepness/2);
end
end
