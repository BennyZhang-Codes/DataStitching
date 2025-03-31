function [f] = plot_grad(Grads,dt)
%PLOT_GRAD
% INPUT:
%   Grads  : [nSample, nTerm] data of dynamic field measurement  
%   dt     : [s] interval time between two points 
% OUTPUT:
%   f      : figure object
% EXAMPLE:
% >> load('grad_1mm.mat');
% >> f = plotSkopeGrads(skopeStitched, dt);
% >> saveas(f, '1', 'png')
% >> exportgraphics(f,'f.png','Resolution',300, 'BackgroundColor', '#1F1F1F')
% created by Jinyuan Zhang, 9/29/2024
    figname = 'skope grads';
    position = [100, 300, 900, 600];
    
    color_dfc       = ["#5f4690", "#1d6996", "#38a6a5", "#0f8554", "#73af48", "#edad08", "#e17c05", "#cc503e", "#94346e"];
    legend_dfc      = ["1", "x", "y", "z", "xy", "zy", "3z² - (x² + y² + z²)", "xz", "x² - y²"];
    color_facecolor = "#FFFFFF";
    color_label     = "#010101";
    font_label      = "Times New Roman";

    f = figure('Name', figname, 'Position', position, 'Color', color_facecolor);
    
    [nPoint, nTerm] = size(Grads);
    t = (1:nPoint) * dt * 1e3; % ms
    
    subplot(3,1,1, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;
    plot(t, real(Grads(:, 1)), 'Color', color_dfc(1));
    ylim([-0.02, 0.02]);
    xlabel('Time [ms]'     , 'Color', color_label, 'fontname', font_label);
    ylabel('0th order [mT]', 'Color', color_label, 'fontname', font_label);
    lgd = legend(legend_dfc(1),"TextColor",color_label,"Box","off","FontName",font_label, "NumColumns",1, "Location","northwest");
    lgd.ItemTokenSize = 10 * ones(1,5);

    subplot(3,1,2, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;
    plot(t, Grads(:, 2), 'Color', color_dfc(2));
    plot(t, Grads(:, 3), 'Color', color_dfc(3));
    plot(t, Grads(:, 4), 'Color', color_dfc(4));
    ylim([-60, 60]);
    xlabel('Time [ms]'       , 'Color', color_label, 'fontname', font_label);
    ylabel('1st order [mT/m]', 'Color', color_label, 'fontname', font_label);
    lgd = legend(legend_dfc(2:4),"TextColor",color_label,"Box","off","FontName",font_label, "NumColumns",3, "Location","northwest");
    lgd.ItemTokenSize = 10 * ones(1,5);

    subplot(3,1,3, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;
    ylim([-2, 2]);
    plot(t, Grads(:, 5), 'Color', color_dfc(5))
    plot(t, Grads(:, 6), 'Color', color_dfc(6))
    plot(t, Grads(:, 7), 'Color', color_dfc(7))
    plot(t, Grads(:, 8), 'Color', color_dfc(8))
    plot(t, Grads(:, 9), 'Color', color_dfc(9))
    xlabel('Time [ms]'         , 'Color', color_label, 'fontname', font_label);
    ylabel('2nd order [mT/m^2]', 'Color', color_label, 'fontname', font_label);
    lgd = legend(legend_dfc(5:9),"TextColor",color_label,"Box","off","FontName",font_label, "NumColumns",5, "Location", "northwest");
    lgd.ItemTokenSize = 10 * ones(1,5);
end





