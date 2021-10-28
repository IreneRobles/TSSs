module TSSs

using DataFrames
using Statistics
using DelimitedFiles
using PyPlot
using Images
using Colors
using ImageMagick
using FileIO
using ColorTypes
using Distances
using ProgressMeter
using ImageFiltering

include("FilesAndImages.jl")
include("ReadOutlines.jl")
include("GeneralFunctions.jl")
include("TSSsingleChannel.jl")
include("Plot.jl")
include("linkTSS.jl")



function quant(tss, tssname,  cell, key)
    try
        tss[cell][tssname][key][1]
        catch; 0
    end
end

function add_tss_quantification(t2,t3,t4,freq_tb; quantmethod = :Int_Int_radious2)
    cells = freq_tb[:Cell]
    
    freq_tb[:locus1_TSS2_size] = [quant(t2, freq_tb[ii, :locus1_tss2], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[:locus2_TSS2_size] = [quant(t2, freq_tb[ii, :locus2_tss2], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[:locus1_TSS3_size] = [quant(t3, freq_tb[ii, :locus1_tss3], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[:locus2_TSS3_size] = [quant(t3, freq_tb[ii, :locus2_tss3], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[:locus1_TSS4_size] = [quant(t4, freq_tb[ii, :locus1_tss4], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[:locus2_TSS4_size] = [quant(t4, freq_tb[ii, :locus2_tss4], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]

    freq_tb
       
end





end # module
