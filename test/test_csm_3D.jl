using Test
using HighOrderMRI
import HighOrderMRI: load_csm
nX              = 200   ;  
nY              = 200   ;
nZ              = 100  ;
nCoil           = 18    ;
nRow            = 3     ;
nCol            = 6     ;
nBlock          = 3     ;
overlap         = 1     ;
relative_radius = 1.5   ;
verbose         = true  ;

@testset "csm 3D" begin
    @testset "csm_Real_32cha 3D" begin
        csm = load_csm(:real_32cha, nX, nY, nZ, nCoil; verbose=verbose);
        @test size(csm) == (nX, nY, nZ, 32)
        @test typeof(csm) <: AbstractArray{<:Number,4}
        # plt_images(abs.(csm[:,:,35,:]); dim=3, nRow=4, nCol=8, title="csm_Real_32cha 3D")
    end
end

