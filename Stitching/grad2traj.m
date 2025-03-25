function [y] = grad2traj(x,dt)
%GRAD2TRAJ Trapezoidal cumulative integration over time
% INPUT:
%   x      : [nSample, nTerm] data   
%   dt     : [s] interval time between two points 
% OUTPUT:
%   y      : cumulative integration over time
% created by Jinyuan Zhang, 22/10/2024
    [nSample, nTerm] = size(x);
    y =  (x(2:end, :) + x(1:end-1, :)) .* (dt / 2);
    y = cumsum(y, 1);
    y = [zeros(1,nTerm); y];
end

