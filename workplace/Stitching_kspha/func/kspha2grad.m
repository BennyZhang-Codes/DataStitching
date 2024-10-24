function [dt, nSeg, nSampleAllSeg, bfield_stitched] = kspha2grad(params, delay_offset)
%kspha2grad Load skope recording and convert to gradient from bfield
% INPUT:
%   params          : struct
%   delay_offset    : delay_offset * dt [s], fine tuning of the delay time manually.
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
        fprintf("Stitching method\n");
        nSampleAllSeg = round((trigDelay(2:end)-trigDelay(1:end-1))/dt);
        nSampleAllSeg = [nSampleAllSeg; round((readoutGradDur-trigDelay(end))/dt)];
    else
        fprintf("Conventional method\n");
        nSampleAllSeg = nSamplePerSeg;
    end
    nSampleStitched = sum(nSampleAllSeg);
    %% skope_phaseCoeffs2Bfield
    clear bfield
    delayshift = gradFreeDelay- scan.k.extTrigDelay + delay_offset*dt - 0.5*dt;  % 0.5*dt for the shift caused by "diff" in "deriveBfieldFromPhase" function
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

    %% Stitching of all the segments
    bfield_stitched = zeros(nSampleStitched,nTerm);
    for seg = 1:nSeg
        idx_e = sum(nSampleAllSeg(1:seg));
        idx_s = sum(nSampleAllSeg(1:seg)) - nSampleAllSeg(seg) + 1; % fprintf("%d : %d\n", idx_s, idx_e);
        bfield_stitched(idx_s:idx_e,:) = bfield(1:nSampleAllSeg(seg),:,1,seg);
    end
    %% 
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