function [traj_interp] = interp1_TrajTime1(traj, dt, delTime, datatime)
% Interpolates trajectory along the 1st dimension 
% INPUT:
%   dt        : time delay between input samples
%   delTime   : a time delay that is effectively added to datatime (or subtracted from the times of the input signals)
%   datatime  : Output time points
% OUTPUT:
%    Output is at samples defined by datatime
%
%   Time for trajectory input is presumed to start at 0.
%
% created by Jinyuan Zhang 2024
    trajTim =  - delTime + dt * (0:size(traj,1)-1); 
    traj_interp = interp1(trajTim, traj, datatime, 'makima',0); 
end
