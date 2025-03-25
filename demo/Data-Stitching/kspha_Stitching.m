close all; clear;

addpath('../../Stitching');
addpath(genpath('PATH_TO_SKOPE_AqSysDataImport'));

root         = '.';
seqsStitched = ["7T_1p0_200_r4_seg4_max51.seq", ];
idsStitched  = [25,];
seqsStandard = ["7T_1p0_200_r4_seg1_max51.seq", ];
idsStandard  = [24];
names        = ["1p0_200_r4", ];

for idx = 1:1
    run(root, seqsStitched(idx), seqsStandard(idx), idsStitched(idx), idsStandard(idx), names(idx))
end


function run(root, seqStitched, seqStandard, idStitched, idStandard, name)
    color_facecolor = "#FFFFFF";
    paramsStitched = struct();
    paramsStitched.folder       = fullfile(root, 'rawdata_skope');
    paramsStitched.id           = idStitched;
    paramsStitched.seqfile      = fullfile(root, 'seq', seqStitched);
    paramsStitched.delay_offset = 0*1e-6;
    paramsStitched.mode         = 'variable';  % variable
    
    paramsStandard = struct();
    paramsStandard.folder       = fullfile(root, 'rawdata_skope');
    paramsStandard.id           = idStandard;
    paramsStandard.seqfile      = fullfile(root, 'seq', seqStandard);
    paramsStandard.delay_offset = 0*1e-6;
    paramsStandard.mode         = 'variable';  % variable


    %% Stitching measurement
    [dt, nSegStitched, nSampleAllSegStitched, bfieldStitched] = kspha2grad(paramsStitched);
    fBfieldStitched = plot_grad(bfieldStitched, dt); 
    exportgraphics(fBfieldStitched, fullfile(root, 'result', sprintf('7T_%s_Stitched_Bfield.png', name)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
    
    [delayStitched, ~, ~, ~, ksphaStitched] = ksphaStitch(paramsStitched);
    fKsphaStitched  = plot_kspha(ksphaStitched, dt); 
    exportgraphics(fKsphaStitched,  fullfile(root, 'result', sprintf('7T_%s_Stitched_kspha.png', name)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
    
    
    %% Conventional measurement
    [dt, nSegStandard, nSampleAllSegStandard, bfieldStandard] = kspha2grad(paramsStandard);
    fBfieldStandard = plot_grad(bfieldStandard, dt); 
    exportgraphics(fBfieldStandard, fullfile(root, 'result', sprintf('7T_%s_Standard_Bfield.png', name)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
    [delayStandard, ~, ~, ~, ksphaStandard] = ksphaStitch(paramsStandard);
    fKsphaStandard  = plot_kspha(ksphaStandard, dt); 
    exportgraphics(fKsphaStandard,  fullfile(root, 'result', sprintf('7T_%s_Standard_kspha.png', name)), 'Resolution', 300, 'BackgroundColor', color_facecolor);
    
    fksphadiff = plot_kspha_xy_comparison({ksphaStandard, ksphaStitched}, ["Standard", "Stitched"]);
    exportgraphics(fksphadiff,  fullfile(root, 'result', sprintf('7T_%s_ksphadiff.png', name)), 'Resolution', 300, 'BackgroundColor', color_facecolor);

    %% Saving
    % save(fullfile(root, 'dfc_2nd', sprintf('7T_%s.mat', name)), 'delayStitched', 'nSampleAllSegStitched', 'bfieldStitched', 'ksphaStitched', 'dt');
    % save(fullfile(root, 'dfc_2nd', sprintf('7T_%s.mat', name)), 'delayStandard', 'nSampleAllSegStandard', 'bfieldStandard', 'ksphaStandard', '-append');
    % close all;
end


