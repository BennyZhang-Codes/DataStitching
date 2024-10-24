close all;
clear;

root = 'E:\skope\AqSysDataImport\';addpath(root);skopepath = fullfile(root,'24_09_24_spiral_isocenter');
r = 4;
idStitched = 8;  % 4   8-4segs
seqStitched = sprintf('xw_sp2d_7T-1mm-200-r%d-4segs_woRF.seq', r);
seqfileStitched = fullfile(skopepath, seqStitched);
plotfileStitched = fullfile(skopepath, sprintf('r%d_Stitched.png', r));
%% segmented measurement
params = struct();
params.folder = skopepath;
params.id = idStitched;
params.seqfile = seqfileStitched;
params.mode = 'variable';  % variable


% [gradStitched, dt, nSegStitched, gradSegStitched, ntStitched, skopeStitched] = skope2grad1(params);
% fStitched = plotSkopeGrads(skopeStitched, dt); exportgraphics(fStitched, plotfileStitched,'Resolution',300, 'BackgroundColor', '#1F1F1F')

folder = params.folder;
id = params.id;
seqName = params.seqfile;
mode = params.mode;

% 1. read *.seq, get some corresponding parameters
myseq = mr.Sequence();
myseq.read(seqName);
gradFreeDelay  = myseq.getDefinition('skope_gradFreeDelay');
nSeg           = myseq.getDefinition('skope_nSegments2measure');
gradRasterTime = myseq.getDefinition('GradientRasterTime');
nPrescan       = myseq.getDefinition('skope_nPrescans');
trigDelay      = myseq.getDefinition('skope_triggerDelays');
readoutGradDur = myseq.getDefinition('readoutGradientDuration');
segDur = myseq.getDefinition('skope_segmentDuration');
if isempty(segDur)
    % mode = 'variable';
    segDur = myseq.getDefinition('skope_maxSegmentDuration');
end

% 2. read 
scan = AqSysData(folder, id);
dt = scan.k.tDwell;

% [nSample,nChannel,nInterleave,nDynamic]
kspha = scan.getData('kspha', [], [], [], nPrescan+1:nPrescan+nSeg); 
[nSample,nTerm,nInterleave,nDynamic] = size(kspha);
% figure, plot(kspha(:,1:4)); xlabel('Samples'), ylabel('k_0 [rad], k_x_y_z [rad/m]'), legend('k_0','k_x','k_y','k_z'), title('k vs time')
% figure, plot(kspha(:,2),kspha(:,3)); xlabel('kx [rad/m]'), ylabel('ky [rad/m]'), axis equal, title('k parametric')

tshift = round((gradFreeDelay- scan.k.extTrigDelay)./dt);
nSamplePerSeg = round(segDur/dt); % number of samples per segmentation.


%read raw data [nSample,nChannel,nInterleave,nDynamic]
raw = scan.getData('raw', [], [], [], nPrescan+1:nPrescan+nSeg);



ProbeData.B0 = scan.B0;
ProbeData.dataRaw = raw;
ProbeData.dt = scan.k.tDwell;
ProbeData.fieldOffsets = scan.probeFieldOffset;
ProbeData.gammaMRI = scan.gammas.mrSystem;
ProbeData.gammaProbes = scan.gammas.probes;
ProbeData.probePositions = scan.probePositions;

% Perform a ordinary least squares solution. 
opt.fitOrder = 2;          % Fit to 2nd order spherical harmonics
opt.fitDoSteps = false;    % Stepwise correction
opt.fitDoWeights = false;  % Weighted least square correction (based on probe distance from isocenter)
opt.fitDoResAdj = false;   % Weighted residual correction
[phs_spha_ls, phs_conc_ls] = harmonicsFromRaw(ProbeData.probePositions,ProbeData.dataRaw,ProbeData.dt,ProbeData.fieldOffsets,ProbeData.gammaProbes,ProbeData.gammaMRI,ProbeData.B0,opt);



%% Plotting
figure;
NP = size(phs_spha_ls,2);
nr = floor(sqrt(NP));
nc = ceil(NP/nr);
t = dt*(1:size(phs_spha_ls,1))*1000; % [ms]
ylabelsAll = {'rad', 'rad/m', 'rad/m^2', 'rad/m^3'};
for np = 1:NP
    subplot(nr,nc,np);
    plot(t,phs_spha_ls(:,np, 1));
    hold('all')
    plot(t,kspha(:,np,1));
    xlabel('time (ms)')
    if np > 9
        ylabel(ylabelsAll{4});
    elseif np > 4
        ylabel(ylabelsAll{3});
    elseif np > 1
        ylabel(ylabelsAll{2});
    else 
        ylabel(ylabelsAll{1});
    end
end
legend('least squares', 'kspha')