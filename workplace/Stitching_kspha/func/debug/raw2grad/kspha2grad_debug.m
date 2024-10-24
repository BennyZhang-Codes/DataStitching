close all;
clear;

root = 'E:\skope\AqSysDataImport\';addpath(root);skopepath = fullfile(root,'24_09_24_spiral_isocenter');
r = 4;
idStitched =8;  % 4   8-4segs
seqStitched = sprintf('xw_sp2d_7T-1mm-200-r%d-4segs_woRF.seq', r);
seqfileStitched = fullfile(skopepath, seqStitched);
plotfileStitched = fullfile(skopepath, sprintf('r%d_Stitched.png', r));
%% segmented measurement
params = struct();
params.folder = skopepath;
params.id = idStitched;
params.seqfile = seqfileStitched;

folder = params.folder;     % root directory
id = params.id;             % ID of the skope recording
seqName = params.seqfile;   % full path of *.seq file

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

% 2. read skope raw data
scan = AqSysData(folder, id);
dt = scan.k.tDwell;
nSamplePerSeg = round(segDur/dt); % number of samples per segmentation.


% [nSample,nChannel,nInterleave,nDynamic]
kspha = scan.getData('kspha', [], [], [], nPrescan+1:nPrescan+nSeg); 
[nSample,nTerm,nInterleave,nDynamic] = size(kspha);
%% determination of nSampleAllSeg and nSampleStitched
if length(trigDelay) > 1   % if length(trigDelay) > 1, stitching method, otherwise, conventional method
    nSampleAllSeg = round((trigDelay(2:end)-trigDelay(1:end-1))/dt);
    nSampleAllSeg = [nSampleAllSeg; round((readoutGradDur-trigDelay(end))/dt)];
else
    nSampleAllSeg = nSamplePerSeg;
end
nSampleStitched = sum(nSampleAllSeg);
%% skope_phaseCoeffs2Bfield
clear bfield;

delayshift = gradFreeDelay- scan.k.extTrigDelay - 0.5*dt;
bfield = zeros(nSample-1,nTerm,nInterleave,nDynamic);
for dynamic = 1: nDynamic 
    for interleave =1 : nInterleave
        for term = 1: nTerm 
            bfield(:,term,interleave,dynamic) = deriveBfieldFromPhase(kspha(:,term,interleave,dynamic),dt,gradRasterTime,0);
            bfield(:,term,interleave,dynamic) = interp1_TrajTime(bfield(:,term,interleave,dynamic), dt, delayshift);
        end
    end
end
bfield = bfield(1:nSamplePerSeg,:,:,:); % [nSamplePerSeg,nTerm,nInterleave,nDynamic]
bfield = mr.convert(bfield,'Hz/m','mT/m','gamma',scan.gammas.mrSystem);

bfield_stitched = zeros(nSampleStitched,nTerm);
for seg = 1:nSeg
    idx_e = sum(nSampleAllSeg(1:seg));
    idx_s = sum(nSampleAllSeg(1:seg)) - nSampleAllSeg(seg) + 1; % fprintf("%d : %d\n", idx_s, idx_e);
    bfield_stitched(idx_s:idx_e,:) = bfield(1:nSampleAllSeg(seg),:,1,seg);
end
plot_grad(bfield_stitched, dt);
%%
bfield_traj = grad2traj(bfield_stitched(1:end,:), dt) * 1e-3 * mr.opts().gamma*2*pi;
kspha_traj = interp1_TrajTime(kspha(:, :,1,1), dt, gradFreeDelay- scan.k.extTrigDelay);
figure,subplot(1,2,1),hold on;
plot(kspha_traj(1:nSamplePerSeg,2));
plot(bfield_traj(1:nSamplePerSeg, 2));
subplot(1,2,2)
plot(kspha_traj(1:nSamplePerSeg,2)-bfield_traj(1:nSamplePerSeg, 2));


