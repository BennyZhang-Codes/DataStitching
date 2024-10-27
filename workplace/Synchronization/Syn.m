close all;
clear;

addpath(genpath('/home/jyzhang/Desktop/synchronization/matmri'))
root = '/home/jyzhang/Desktop/synchronization/Stitching';
addpath(fullfile(root,'func'))
addpath(fullfile(root,'Synchronization'))

sub_folder = '/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo';
cd(sub_folder);
%%
% ΔB0, Coil-Sensitivity, Mask
MID = 66;
grefile = dir(fullfile("syn", sprintf('syn*MID*%d*FID*gres6.mat', MID))).name;
load(fullfile("syn", grefile));

% MRI signal data: 
r = 1;
MRIfile = dir(fullfile("syn", sprintf('syn*MID*FID*r%d*.mat', r))).name;
load(fullfile("syn", MRIfile), 'data', 'FOV', 'matrixSize');

% Dynamic Field  Changes: kspha, dt
DFCfile = dir(fullfile("dfc", sprintf('xw_sp2d_7T*r%d*.mat', r))).name;
load(fullfile("dfc", DFCfile), 'ksphaStitched', 'ksphaStandard', 'dt', 'delayStitched', 'delayStandard');

%% preparation of params
params.b0map         = -2 * pi * b0;  % [Hz] -> [rad/s] figure,subplot(1,1,1); imagesc(b0); axis('image'); axis('off'); colormap('gray');
params.csm           = csm;
params.mask          = mask;

params.data          = data;
params.FOV           = FOV;
params.matrixSize    = matrixSize;

params.dt            = dt;
params.kspha         = ksphaStandard;
params.dataStartTime = delayStandard-1*dt;

%% do Synchronization
[ksphaStandard_syn, syndelStandard, ~,~] = Syn_kspha(params);
save(fullfile("syn", MRIfile), "ksphaStandard_syn", "syndelStandard", "-append")

params.kspha         = ksphaStitched;
params.dataStartTime = delayStitched-1*dt;
[ksphaStitched_syn, syndelStitched, ~,~] = Syn_kspha(params);
save(fullfile("syn", MRIfile), "ksphaStitched_syn", "syndelStitched", "-append")

