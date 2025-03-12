nSample = 100
nTerm   = 9
dt      = 1.    # dwell time for trajectory is 1
T       = Float64


datatime = T.(collect(0:2.5:nSample-1)); # dwell time for data is 2.5
kspha  = T.(sin.((collect(1:nSample)) ./ nSample .* π .* collect(1:9)'));

kspha_del     = InterpTrajTime(kspha    , dt,  1*dt);
kspha_del_inv = InterpTrajTime(kspha_del, dt, -1*dt); 
# plt_ksphas([kspha, kspha_del_inv], dt)

@test kspha[2:end-2,:] ≈ kspha_del_inv[2:end-2,:]