using Test
using HighOrderMRI
import HighOrderMRI: load_csm
nX              = 200   ;  
nY              = 200   ;
nCoil           = 18    ;
nRow            = 3     ;
nCol            = 6     ;
nBlock          = 3     ;
overlap         = 1     ;
relative_radius = 1.5   ;
verbose         = false ;

@testset "csm" begin
    @testset "csm_Fan_binary" begin
        csm = load_csm(:fan, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Real,3}
        # plt_images(csm; dim=3, nRow=nRow, nCol=nCol, title="csm_Fan_binary")
    end
    @testset "csm_Rect_binary" begin
        csm = load_csm(:rect, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Real,3}
        # plt_images(csm; dim=3, nRow=nRow, nCol=nCol, title="csm_Rect_binary")
    end
    @testset "csm_Rect_gaussian" begin
        csm = load_csm(:rect_gaussian, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(  abs.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Rect_gaussian_mag")
        # plt_images(angle.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Rect_gaussian_pha")
    end
    @testset "csm_Birdcage" begin
        csm = load_csm(:birdcage, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(abs.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Birdcage")
    end
    @testset "csm_Real_32cha" begin
        csm = load_csm(:real_32cha, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, 32)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(abs.(csm); dim=3, nRow=4, nCol=8, title="csm_Real_32cha")
    end
    @testset "csm_Gaussian_grid" begin
        csm = load_csm(:gaussian_grid, nX, nY, nCoil; nRow=nRow, nCol=nCol, overlap=overlap, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(abs.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Gaussian_grid")
    end
    @testset "csm_Gaussian_grid_block" begin
        csm = load_csm(:gaussian_grid_block, nX, nY, nCoil; nRow=nRow, nCol=nCol, nBlock=nBlock, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(abs.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Gaussian_grid_block")
    end
    @testset "csm_Gaussian_grid_block_pha" begin
        csm = load_csm(:gaussian_grid_block_pha, nX, nY, nCoil; use_gpu=false, nRow=nRow, nCol=nCol, nBlock=nBlock, relative_radius=relative_radius, verbose=verbose)
        @test size(csm) == (nX, nY, nCoil)
        @test typeof(csm) <: AbstractArray{<:Number,3}
        # plt_images(  abs.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Gaussian_grid_block_pha mag")
        # plt_images(angle.(csm); dim=3, nRow=nRow, nCol=nCol, title="csm_Gaussian_grid_block_pha pha")
    end
end

