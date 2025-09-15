using HighOrderMRI
using PyPlot; PyPlot.pygui(true)
using MAT

path = "demo3D";

seq = read_seq("$(path)/Spiral_2p0x2p0x2p0_100x100x35_1x1_int1_9p5ms_tr500_te5p0.seq")
seq.GR[1,:] = -seq.GR[1,:]

hoseq = HO_Sequence(seq);

params = matread("$(path)/Spiral_2p0x2p0x2p0_100x100x35_1x1_int1_9p5ms_tr500_te5p0.mat")
dt_adc         = params["skope"]["dt_adc"];
nSample_adc    = Int64(params["skope"]["nSample_adc"]);
dur_adc        = nSample_adc*dt_adc;
gradRasterTime = 10e-6;
nGrad_adc      = Int64(round(dur_adc/gradRasterTime));

skope = matread("$(path)/Spiral_2p0x2p0x2p0_100x100x35_1x1_int1_9p5ms_tr500_te5p0_skope.mat")
g_skope = skope["g_skope"] * 1e-3; # in T/m
t_g     = vec(skope["t_g"]);


nBlock = size(seq.GR, 2)
iadc = 1;
for iblock = 1:nBlock
    if seq.GR[1, iblock].delay > 0 && is_ADC_on(seq[iblock])
        g = InterpTrajTime(g_skope[iadc, :, :], gradRasterTime, -0.9*gradRasterTime, t_g)

        nGrad_Gz_pre  = 34;
        nGrad_Spoiler = 300;
        nGrad_adc     = Int64(round(dur_adc/gradRasterTime));
        nGrad_delay   = Int64(round(seq.GR[1, iblock].delay/gradRasterTime));
    
        for i = 1:16
            if i == 4
                hoseq.GR_dfc[i, iblock-1].T    = (nGrad_Gz_pre-1)*gradRasterTime;
                hoseq.GR_dfc[i, iblock-1].A    = g[1:nGrad_Gz_pre, i];
                hoseq.GR_dfc[i, iblock-1].rise = gradRasterTime/2;
                hoseq.GR_dfc[i, iblock-1].fall = gradRasterTime/2;
            end

            hoseq.GR_dfc[i, iblock].T     = (nGrad_adc-1)*gradRasterTime;
            hoseq.GR_dfc[i, iblock].delay = nGrad_delay*gradRasterTime;
            hoseq.GR_dfc[i, iblock].A     = g[nGrad_Gz_pre+nGrad_delay+1:end, i];
            hoseq.GR_dfc[i, iblock].rise  = gradRasterTime/2;
            hoseq.GR_dfc[i, iblock].fall  = gradRasterTime/2;
            if i == 2 || i == 3
                hoseq.GR_dfc[i, iblock].T = (nGrad_adc+nGrad_Spoiler-1)*gradRasterTime;
                hoseq.GR_dfc[i, iblock].A = [g[nGrad_Gz_pre+nGrad_delay+1:end, i]; seq.GR[i-1, iblock].A[end-nGrad_Spoiler+1:end]]
            end
        end
        iadc += 1
    else
        hoseq.GR_dfc[2:4, iblock]  = deepcopy(hoseq.SEQ.GR[1:3, iblock]); 
    end
end
plt_seq(hoseq)
plot_seq(hoseq.SEQ)
