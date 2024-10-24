function [dt, nSeg, nSampleAllSeg, bfield_stitched] = kspha2grad0(params)
%kspha2grad0 Load skope recording and convert to gradient from bfield
% INPUT:
%   params          : struct
% OUTPUT:
%   dt              : sampling interval
%   nSeg            : number of segments
%   nSampleAllSeg   : number of samples in each segment
%   bfield_stitched : number of timepoints in every segments
% created by Zihao Zhang, 10/20/2023
% modified by Jinyuan Zhang, 10/19/2024
    folder = params.folder;     % root directory
    id = params.id;             % ID of the skope recording
    seqName = params.seqfile;   % full path of *.seq file
    mode = params.mode;         % 'variable', 'stitched'

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

    % [nSample,nChannel,nInterleave,nDynamic]
    kspha = scan.getData('kspha', [], [], [], [nPrescan+1:nPrescan+nSeg]); 
    [nSample,nTerm,nInterleave,nDynamic] = size(kspha);


    tshift = round((gradFreeDelay- scan.k.extTrigDelay)./dt);
    fprintf('extTrigDelay=%0.2fus\n', 1e6*scan.k.extTrigDelay);
    nSamplePerSeg = round(segDur/dt); % number of samples per segmentation.

    %% skope_phaseCoeffs2Bfield
    for dynamic = 1: nDynamic 
        for interleave =1 : nInterleave
            for channel = 1: nTerm 
                bfield(:,channel,interleave,dynamic) = deriveBfieldFromPhase(kspha(:,channel,interleave,dynamic),dt,gradRasterTime,tshift);
            end
        end
    end
    bfield = bfield(1:nSamplePerSeg,:,:,:); % [nSamplePerSeg,nTerm,nInterleave,nDynamic]

    %% 
    nSampleAllSeg = nSamplePerSeg;
    bfieldSeg = permute(bfield,[1 4 2 3]);  % [nSamplePerSeg,nDynamic,nTerm,nInterleave], nDynamic <=> nSeg
    bfieldSeg = mr.convert(bfieldSeg,'Hz/m','mT/m','gamma',scan.gammas.mrSystem);

    if strcmp(mode, 'variable')
        nSampleAllSeg = round((trigDelay(2:end)-trigDelay(1:end-1))/dt);
        nSampleAllSeg = [nSampleAllSeg; round((readoutGradDur-trigDelay(end))/dt)];
    
        bfield_stitched = [];
        for seg = 1:nSeg
            bfield_stitched = cat(1, bfield_stitched, squeeze(bfieldSeg(1:nSampleAllSeg(seg),seg,:)));
        end
    elseif strcmp(mode, 'stitched')
        bfield_stitched = reshape(bfieldSeg,[],nTerm); % [nSamplePerSeg×nDynamic,nTerm,nInterleave]
    end
    bfield_stitched = bfield_stitched(:,:,1); % [nSamplePerSeg×nDynamic,nTerm,nInterleave] => % [nSamplePerSeg×nDynamic,nTerm]

    fprintf('extTrigDelay   = %0.2f us\n', 1e6*scan.k.extTrigDelay);
    fprintf('dt             = %0.2f us\n', 1e6*dt                 );
    fprintf('nPrescan       = %d\n',       nPrescan               );
    fprintf('nSegment       = %d\n',       nSeg                   );
    fprintf('nInterleave    = %d\n',       nInterleave            );
    fprintf('nDynamic       = %d\n',       nDynamic               );
    fprintf('nSample        = %d\n',       nSample                );
    fprintf('nSamplePerSeg  = %d\n',       nSamplePerSeg          );
    fprintf('triggerDelays  = %0.2f ms\n', 1e3*trigDelay          );
    fprintf('readoutGradDur = %0.2f ms\n', 1e3*readoutGradDur     );
end