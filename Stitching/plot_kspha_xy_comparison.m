function [f] = plot_kspha_xy_comparison(ksphas, labels)
%PLOT_KSPHA_XY
% INPUT:
%   kspha  : [nSample, nTerm, nSeg] kspha data of dynamic field measurement  
%   dt     : [s] interval time between two points 
% OUTPUT:
%   f      : figure object
% EXAMPLE:
% created by Jinyuan Zhang, 22/10/2024
    figname = 'skope kspha xy';
    position = [100, 500, 500, 500];
    
    color_dfc       = ["#5f4690", "#1d6996", "#38a6a5", "#0f8554", "#73af48", "#edad08", "#e17c05", "#cc503e", "#94346e"];
    legend_dfc      = ["1", "x", "y", "z", "xy", "zy", "3z² - (x² + y² + z²)", "xz", "x² - y²"];
    color_facecolor = "#FFFFFF";
    color_label     = "#010101";
    font_label      = "Times New Roman";
    
    f = figure('Name', figname, 'Position', position, 'Color', color_facecolor);
    
    [nPoint, nTerm] = size(ksphas{1});
    
    subplot(1,1,1, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;
    
    for i = 1:length(ksphas)
        plot(ksphas{i}(:, 2), ksphas{i}(:, 3));
    end
    lgd = legend(labels,"TextColor",color_label,"Box","off","FontName",font_label, "NumColumns",1, "Location", "northeast");
    lgd.ItemTokenSize = 10 * ones(1,5);
    xlabel('X [rad/m]', 'Color', color_label, 'fontname', font_label);
    ylabel('Y [rad/m]', 'Color', color_label, 'fontname', font_label);
    axis equal;
end





