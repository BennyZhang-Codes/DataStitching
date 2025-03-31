function [delay, dt, nSeg, nSampleAllSeg, kspha_fg] = ksphaStitch(params)
%ksphaStitch2 Load skope recording and stitch the kspha
% INPUT:
%   params          : struct
%       folder          : root directory
%       id              : ID of the skope recording
%       seqName         : full path of *.seq file
%       delay_offset    : [s], fine tuning of the delay time manually.
% OUTPUT:
%   delay           : delay of the [gradFreeDelay - extTrigDelay + delay_offset], or the start time point in the trajectory
%   dt              : sampling interval
%   nSeg            : number of segments
%   nSampleAllSeg   : number of samples in each segment
%   kspha_Stitched  : number of timepoints in every segments
% created by Jinyuan Zhang, 27/10/2024
    folder = params.folder;     % root directory
    id = params.id;             % ID of the skope recording
    seqName = params.seqfile;   % full path of *.seq file
    delay_offset = params.delay_offset; % [s], fine tuning of the delay time manually.

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
    %% kspha align and stitch
    clear kspha_Stitched bfield_stitched
    delay = gradFreeDelay - scan.k.extTrigDelay + delay_offset;
    delay_rounded = round((delay)/dt);
    
    % bfield_Stitched = [];
    bfield_filter_Stitched = [];
    if nDynamic == 1   % for conventional method
        for interleave = 1 : nInterleave
            k = kspha(:,:,interleave,1);
            k_f = filterPhaseData(kspha(:,:,interleave,1), dt, gradRasterTime);
        end
        % bfield_Stitched = cat(1, bfield_Stitched, diff(k)./ dt);
        bfield_filter_Stitched = cat(1, bfield_filter_Stitched, diff(k_f)./ dt);
    else
        for dynamic = 1: nDynamic 
            clc; fprintf('\rProgress: %d/%d', dynamic, nDynamic);
    
            for interleave = 1 : nInterleave
                % k_f   = kspha(:,:,interleave,dynamic);
                % k_f   = diff(k_f)./ dt;

                k_f = filterPhaseData(kspha(:,:,interleave,dynamic), dt, gradRasterTime);
                k_f = diff(k_f)./ dt;
                if dynamic == 1
                    % k   = k(1:nSampleAllSeg(dynamic)+delay_rounded,:);
                    k_f = k_f(1:nSampleAllSeg(dynamic)+delay_rounded, :);
                elseif dynamic == nDynamic
                    % k   = k(delay_rounded:end,:);
                    k_f = k_f(1+delay_rounded:end, :);
                else
                    % k   = k(delay_rounded:nSampleAllSeg(dynamic)+delay_rounded,:);
                    k_f = k_f(1+delay_rounded:nSampleAllSeg(dynamic)+delay_rounded,:);

                end
                % disp(size(k_f))
                % bfield_Stitched = cat(1, bfield_Stitched, k);
                bfield_filter_Stitched = cat(1, bfield_filter_Stitched, k_f);
            end
        end
    end
    % kspha_g = grad2traj(bfield_Stitched, dt);
    % kspha_g = interp1_TrajTime1(kspha_g, dt, -0.5*dt, dt * (0:size(kspha_g,1)-1));
    % kspha_gf = filterPhaseData(kspha_g, dt, gradRasterTime);
    % kspha_gf = kspha_gf(1:nSampleStitched+delay_rounded, :);
    
    kspha_fg = grad2traj(bfield_filter_Stitched, dt);
    kspha_fg = interp1_TrajTime1(kspha_fg, dt, -0.5*dt, dt * (0:size(kspha_fg,1)-1));
    kspha_fg = kspha_fg(1:nSampleStitched+delay_rounded, :);

    %% 
    fprintf('extTrigDelay   = %0.3f us\n', 1e6*scan.k.extTrigDelay);
    fprintf('delay          = %0.3f us\n', 1e6*delay              );
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