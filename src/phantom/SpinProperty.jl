
function SpinProperty_1p5T(class::Array)
    # Define spin property vectors
    T2 = (class.==23 )*329 .+ #CSF
         (class.==46 )*83  .+ #GM
         (class.==70 )*70  .+ #WM
         (class.==93 )*70  .+ #FAT1
         (class.==116)*47  .+ #MUSCLE
         (class.==139)*329 .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*70  .+ #FAT2
         (class.==232)*329 .+ #DURA
         (class.==255)*70     #MARROW
    T2s =(class.==23 )*58  .+ #CSF
         (class.==46 )*69  .+ #GM
         (class.==70 )*61  .+ #WM
         (class.==93 )*58  .+ #FAT1
         (class.==116)*30  .+ #MUSCLE
         (class.==139)*58  .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*61  .+ #FAT2
         (class.==232)*58  .+ #DURA
         (class.==255)*61  .+ #MARROW
         (class.==255)*70     #MARROW
    T1 = (class.==23 )*2569 .+ #CSF
         (class.==46 )*833 .+ #GM
         (class.==70 )*500 .+ #WM
         (class.==93 )*350 .+ #FAT1
         (class.==116)*900 .+ #MUSCLE
         (class.==139)*569 .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*500 .+ #FAT2
         (class.==232)*2569 .+ #DURA
         (class.==255)*500    #MARROW
    ρ =  (class.==23 )*1   .+ #CSF
         (class.==46 )*.86 .+ #GM
         (class.==70 )*.77 .+ #WM
         (class.==93 )*1   .+ #FAT1
         (class.==116)*1   .+ #MUSCLE
         (class.==139)*.7  .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*.77 .+ #FAT2
         (class.==232)*1   .+ #DURA
         (class.==255)*.77    #MARROW

    T1  = T1  * 1e-3
    T2  = T2  * 1e-3
    T2s = T2s * 1e-3
    return T1, T2, T2s, ρ
end