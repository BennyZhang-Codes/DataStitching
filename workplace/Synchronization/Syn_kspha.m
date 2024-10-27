function [kspha_syn, delFound, im, imFoundDel] = Syn_kspha(params)
%SYN_KSPHA Synchronization between measurement of dynamic field and the MRI
%data
% INPUT:
%   params            : struct
%       b0map         : Delta B0 map
%       csm           : Coil-Sensitivity map
%       mask          : mask 
%       data          : MRI Signal data
%       FOV           : Field-Of-View
%       matrixSize    : matrix size
%       dt            : dwell time of the trajectory
%       kspha         : 
%       dataStartTime : start time of the MRI acquisition relative to the start of trajectory
% OUTPUT:
%   kspha_syn  : 
%   delFound   : 
    b0map = params.b0map;      
    csm   = params.csm;        
    mask  = params.mask;       
    
    data = params.data;       
    FOV = params.FOV;      
    matrixSize = params.matrixSize; 
    
    dt = params.dt;         
    kspha = params.kspha;      
    dataStartTime = params.dataStartTime; 
    
    
    % MRI signal data, create datatime (signal sampling time) [s]
    [nSample, nCha] = size(data);
    datatime = dt*(0:nSample-1);
    
    % Build receiver operator from Coil-Sensitivity Map
    R = rcvrOp(csm,0);  % figure, imshow3(abs(csm),[],[4,8]); title('Coil-Sensitivity map'); colormap((gray(256))); colorbar;
    
    % sampling grid [m]
    nX = double(matrixSize(1)); nY = double(matrixSize(2)); nZ = double(matrixSize(3));
    resX = double(FOV(1)/nX); resY = double(FOV(2)/nY);  % [m]
    [grid.x, grid.y, grid.z] = meshgrid((-nX/2+0.5:nX/2-0.5)*resX, -(-nY/2+0.5:nY/2-0.5)*resY, 0);
    
    % trajectory
    kspha = interp1_TrajTime1(kspha, dt, dataStartTime, dt * (0:size(kspha,1)-1)); % shift for delayStitched (known delay)
    kspha = kspha(:, 1:9);        % up to second order
    kconc = kspha(:, 1:4)*0;       % do not consider conc field
    
    %% Find delay automatically
    fprintf('Finding delay...\n')
    del0 = 0; % starting guess 
    maxNit_cgne = 20;
    delJumpFact = 3;
    numCoarseSearch = 0;
    % if gpuDeviceCount > 0
    %     % findDelAuto determines whether to use gpu based on input data.
    %     data = gpuArray(data);
    % end
    [delFound, delSk_perIt] = findDelAuto(del0,data,dt/dt,kspha,kconc,datatime/dt,grid,b0map*dt,csm,maxNit_cgne,delJumpFact,numCoarseSearch);
    fprintf('Found delay = %.3f us\n',delFound);
    %% interpTrajTime
    kspha_del0 = interp1_TrajTime1(kspha, dt, 0, datatime);
    kconc_del0 = interp1_TrajTime1(kconc, dt, 0, datatime);
    kspha_syn = interp1_TrajTime1(kspha, dt, delFound*dt, datatime);
    kconc_syn = interp1_TrajTime1(kconc, dt, delFound*dt, datatime);
    % Plotting
    plot_kspha_comparison(kspha_del0, kspha_syn, dt);
    
    %% Build sampling object based on trajectory with no delay
    fprintf('Performing recon w/o Synchronization...\n')
    S = sampHighOrder(b0map,datatime(:),kspha_del0', kconc_del0',grid, [], true, true);
    
    opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
    maxIt = 20;
    
    tic; [im, resvec, mse, xnorm, xdiff] = cgne(opFunc,data,[],maxIt); toc;
    figure,subplot(1,1,1); imagesc(abs(im),[0 0.006]); title('w/o Synchronization'); axis('image'); axis('off'); colormap('gray');
    
    %% Determine image with optimal delay
    fprintf('Performing recon w/ Synchronization...\n')
    S = sampHighOrder(b0map,datatime(:),kspha_syn', kconc_syn',grid, [], true, true);
    opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
    tic; imFoundDel = cgne(opFunc,data,[],maxIt); toc;
    figure,subplot(1,1,1); imagesc(abs(imFoundDel),[0 0.006]); title('w/ Synchronization'); axis('image'); axis('off'); colormap('gray');
    
    %% Plotting
    figure;
    subplot(1,3,1); imagesc(abs(im)                ,[0 0.006]); title('w/o Synchronization'); axis('image'); axis('off'); colormap('gray');
    subplot(1,3,2); imagesc(abs(imFoundDel)        ,[0 0.006]); title('w/ Synchronization'); axis('image'); axis('off'); colormap('gray');
    subplot(1,3,3); imagesc(abs(im)-abs(imFoundDel),[0 0.0001]); title('difference map'); axis('image'); axis('off'); colormap('gray');
end

