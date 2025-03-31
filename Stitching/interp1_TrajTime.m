function [traj_interp] = interp1_TrajTime(traj, dt, delTime)
% Interpolates trajectory along the 1st dimension 
% INPUT:
%   dt        : time delay between input samples
%   delTime   : a time delay that is effectively added to datatime (or subtracted from the times of the input signals)
%   datatime  :
% OUTPUT:
%    Output is at samples defined by datatime
%
%   Time for trajectory input is presumed to start at 0.
%
% created by Jinyuan Zhang 2024
    datatime = dt * (0:size(traj,1)-1);
    trajTim = datatime - delTime; 
    traj_interp = interp1(trajTim, traj, datatime, 'makima',0); 
end
