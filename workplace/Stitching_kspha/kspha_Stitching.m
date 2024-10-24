close all;
clear;

root = 'E:\skope\AqSysDataImport\';addpath(root);skopepath = fullfile(root,'24_09_24_spiral_isocenter');
r = 4;
nameStitched = dir(fullfile(skopepath, sprintf('*seg4r%d*.raw', r))).name;
nameStandard = dir(fullfile(skopepath, sprintf('*seg1r%d*.raw', r))).name;

paramsStitched = struct();
paramsStitched.folder = skopepath;
paramsStitched.id = str2num(nameStitched(1));
paramsStitched.seqfile = fullfile(skopepath, sprintf('xw_sp2d_7T-1mm-200-r%d-4segs_woRF.seq', r));
paramsStitched.mode = 'variable';  % variable

paramsStandard = struct();
paramsStandard.folder = skopepath;
paramsStandard.id = str2num(nameStandard(1));
paramsStandard.seqfile = fullfile(skopepath, sprintf('xw_sp2d_7T-1mm-200-r%d_woRF.seq', r));
paramsStandard.mode = 'variable';  % variable
%% setting for figure
color_facecolor = "#1F1F1F";
%% Stitching measurement
[dt, nSegStitched, nSampleAllSegStitched, bfieldStitched] = kspha2grad(paramsStitched, (0.5-0.192)*1e-6);
fBfieldStitched = plot_grad(bfieldStitched, dt); exportgraphics(fBfieldStitched, fullfile(skopepath, sprintf('r%d_Bfield_Stitched.png', r)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
% [~, ~, ~, ksphaStitched] = ksphaStitch(paramsStitched, 0);
% fKsphaStitched  = plot_kspha(ksphaStitched, dt); exportgraphics(fKsphaStitched,  fullfile(skopepath, sprintf('r%d_Kspha_Stitched.png',  r)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
% %% Conventional measurement
% [dt, nSegStandard, nSampleAllSegStandard, bfieldStandard] = kspha2grad(paramsStandard, 0);
% fBfieldStandard = plot_grad(bfieldStandard, dt); exportgraphics(fBfieldStandard, fullfile(skopepath, sprintf('r%d_Bfield_Standard.png', r)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
% [~, ~, ~, ksphaStandard] = ksphaStitch(paramsStandard, 0);
% fKsphaStandard  = plot_kspha(ksphaStandard, dt); exportgraphics(fKsphaStitched,  fullfile(skopepath, sprintf('r%d_Kspha_Standard.png',  r)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
% %% Saving
% save(fullfile(skopepath, sprintf('xw_sp2d_7T-1mm-200-r%d.mat', r)), 'dt', 'nSampleAllSegStitched', 'nSampleAllSegStandard', 'bfieldStitched', 'bfieldStandard', 'ksphaStitched', 'ksphaStandard');
% close all;
%% compare the diff
% plot_grad(bfieldStandard-bfieldStitched, dt);
% plot_kspha(ksphaStitched-ksphaStandard,dt)




%%
% kStandard = grad2traj(bfieldStandard*1e-3*mr.opts().gamma*2*pi, dt);
% kStitched = grad2traj(bfieldStitched*1e-3*mr.opts().gamma*2*pi, dt);
% plot_kspha(ksphaStitched-kStitched, dt);
% % plot_kspha(kStitched, dt);
% % plot_kspha(ksphaStandard-ksphaStitched, dt);
% % plot_kspha(ksphaStandard-kStandard, dt);
% % plot_kspha(ksphaStitched-kStitched, dt);
% % plot_kspha(ksphaStandard-kStitched, dt);
% % plot_kspha(ksphaStandard-ksphaStitched, dt);
% 
% 
% % figure, hold on;
% % plot(kStitched(:,2));
% % plot(ksphaStitched(:,2))