

function count_tss_in_image(imagedict)
    sum([length(keys(imagedict[ii])) for ii in keys(imagedict)])
end



function meassure_tss(tss_outline, tss_image; xy = 0.189, zx = 0.5)
    
    tss_outline_n = count_tss_in_image(tss_outline)
    img = []
    
    if tss_outline_n > 0
    #convert image to grayscale
    # Do not load image if there is no TSS
    img = read_tiff_as_gray(tss_image)
    end
    ## 19thApril2021 FALSE THE PROBLEM WAS X and Y where swithched
    ## Center xy is the center of the TSS

    
    for cell in keys(tss_outline)
        for tss in keys(tss_outline[cell])
            
            x = Int.(tss_outline[cell][tss][:x])
            y = Int.(tss_outline[cell][tss][:y])
     
            tss_image = []
            
            try
              tss_image = img[minimum(x)-3:maximum(x)+3,minimum(y)-3:maximum(y)+3,  :]
            catch 
              tss_image = img[minimum(x):maximum(x),minimum(y):maximum(y),  :]
            end
                
             #patch_size = [3, 3, 3]
            #tss_image = mapwindow(median, tss_image, patch_size)

         max_indexes = [ii for ii in Tuple(findmax(tss_image)[2])]

            tss_outline[cell][tss][:Avg_Int] =  [float(mean(tss_image))]
            tss_outline[cell][tss][:Int_Int] = [float(sum(tss_image))]
            tss_outline[cell][tss][:Max_Int] = [float(maximum(tss_image))]
            

            
            radious1 = 0.
            radious2 = 0.
            radious3 = 0.
            
            try 
                rad = 1
                radious1 = TSSs.int_brightest_pixel(tss_image, radious = 1)
                catch; radious1 = 0.; end
         
            try 
                rad = 2;  
                radious2 = TSSs.int_brightest_pixel(tss_image, radious = 2)
                catch; radious2 = 0.; end
            try 
                rad = 3; radious3 = TSSs.int_brightest_pixel(tss_image, radious = 3)
                catch; radious3 = 0.; end
        
            tss_outline[cell][tss][:Int_Int_radious1] = [radious1]
            tss_outline[cell][tss][:Int_Int_radious2] = [radious2]
            tss_outline[cell][tss][:Int_Int_radious3] = [radious3]
            tss_outline[cell][tss][:z] = [max_indexes[3]]
            
            tss_outline[cell][tss][:z_] = [max_indexes[3]*zx]
            tss_outline[cell][tss][:y_] = [sum([y[1],max_indexes[2]])*xy]
            tss_outline[cell][tss][:x_] = [sum([x[1],max_indexes[1]])*xy]
            
            
            
            
            
        end
    end
    
    return tss_outline
    
end

function meassure_tss(tss, imagfolder, pat, channel; xy = 0.189, zx = 0.5)
    image_ts = find_image(imagfolder, pat, channel)
    tss = meassure_tss(tss, image_ts,  xy = xy, zx = zx)
    return tss
end


function analysis_singlechannel(tss; image = "No", tss_name = :TSS, xy_center = true)
    cells = collect(keys(tss))
    n_tss = [length(keys(tss[cell])) for cell in cells]
    bf_image = sum(n_tss) /2/length(cells)
    
    locus1 = []
    locus1_int1 = []
    locus1_int2 = [] 
    locus1_int3 = []
    locus1_x = [] 
    locus1_y = []
    locus1_z = []
    locus2 = []
    locus2_int1 = []
    locus2_int2 = [] 
    locus2_int3 = []
    locus2_x = [] 
    locus2_y = []
    locus2_z = []
    locus_distance = []
    
    x = :x; y = :y; z = :z
    if xy_center
        x = :x_; y = :y_; z = :z_
    end
    
    for cell in cells
        tsss = keys(tss[cell])
    
   if length(tsss) == 0
                push!(locus1, 0)
                push!(locus1_int1, 0)
                push!(locus1_int2, 0)
                push!(locus1_int3, 0)
                push!(locus1_x, 0)
                push!(locus1_y, 0)
                push!(locus1_z, 0)
                push!(locus2, 0)
                push!(locus2_int1, 0)
                push!(locus2_int2, 0)
                push!(locus2_int3, 0)
                push!(locus2_x, 0)
                push!(locus2_y, 0)
                push!(locus2_z, 0)
                push!(locus_distance, 0)
        elseif length(tsss) == 1
                ts = tss[cell][collect(tsss)[1]]
                push!(locus1, collect(tsss)[1])
                push!(locus1_int1, ts[:Int_Int_radious1][1])
                push!(locus1_int2, ts[:Int_Int_radious2][1])
                push!(locus1_int3, ts[:Int_Int_radious3][1])
                push!(locus1_x, ts[x][1])
                push!(locus1_y, ts[y][1])
                push!(locus1_z, ts[z][1])
                push!(locus2, 0)
                push!(locus2_int1, 0)
                push!(locus2_int2, 0)
                push!(locus2_int3, 0)
                push!(locus2_x, 0)
                push!(locus2_y, 0)
                push!(locus2_z, 0)
                push!(locus_distance, 0)
            elseif length(tsss) == 2
                ts = tss[cell][collect(tsss)[1]]
                ts1_x = ts[:x_][1]
                ts1_y = ts[:y_][1]
                ts1_z = ts[:z][1]
                push!(locus1, collect(tsss)[1])
                push!(locus1_int1, ts[:Int_Int_radious1][1])
                push!(locus1_int2, ts[:Int_Int_radious2][1])
                push!(locus1_int3, ts[:Int_Int_radious3][1])
                push!(locus1_x, ts[x][1])
                push!(locus1_y, ts[y][1])
                push!(locus1_z, ts[z][1])
                ts = tss[cell][collect(tsss)[2]]
                ts2_x = ts[x][1]
                ts2_y = ts[y][1]
                ts2_z = ts[z][1]
                push!(locus2, collect(tsss)[2])
                push!(locus2_int1, ts[:Int_Int_radious1][1])
                push!(locus2_int2, ts[:Int_Int_radious2][1])
                push!(locus2_int3, ts[:Int_Int_radious3][1])
                push!(locus2_x, ts[x][1])
                push!(locus2_y, ts[y][1])
                push!(locus2_z, ts[z][1])
                push!(locus_distance, euclidean([ts1_x, ts1_y, ts1_z], [ts2_x, ts2_y, ts2_z]))
        end
                        
                    
    end
    
 

    
    
    new_df = DataFrame(
        Dict(
            :Cell => cells,
            Symbol(string("N_", tss_name)) => n_tss,
            :Image => [image for ii in cells],
            :Image_BF => [bf_image for ii in cells],
            
            Symbol(string("locus1_", tss_name)) =>  locus1, 
            Symbol(string("locus1_x_", tss_name)) =>  locus1_x, 
            Symbol(string("locus1_y_", tss_name)) =>  locus1_y, 
            Symbol(string("locus1_z_", tss_name)) =>  locus1_z, 
            Symbol(string("locus1_int1_", tss_name)) =>  locus1_int1,
            Symbol(string("locus1_int2_", tss_name)) =>  locus1_int2,
            Symbol(string("locus1_int3_", tss_name)) =>  locus1_int3,
            Symbol(string("locus2_", tss_name)) =>  locus2,
            Symbol(string("locus2_x_", tss_name)) =>  locus2_x, 
            Symbol(string("locus2_y_", tss_name)) =>  locus2_y, 
            Symbol(string("locus2_z_", tss_name)) =>  locus2_z, 
            Symbol(string("locus2_int1_", tss_name)) =>  locus2_int1,
            Symbol(string("locus2_int2_", tss_name)) =>  locus2_int2,
            Symbol(string("locus2_int3_", tss_name)) =>  locus2_int3,
            Symbol(string("locus_dist_", tss_name)) =>  locus_distance

        
        )
    )
end

function calculate_bf(t, tname; limit = 0)
    new_df = DataFrames.DataFrame()
    
    t[:Sample] = t[:Sample] 
    samples = unique(t[:Sample])
    new_df[:Sample] = samples
    new_df[:Genotype] = [split(ii, "_")[1] for ii in samples]
    new_df[:Timepoint] = [split(ii, "_")[2] for ii in samples]
    
    function calculate(t, tname)

    new_df[:N_cells] = [sum(t[:Sample] .== ii) for ii in samples]
    sp_sam = [t[t[:Sample] .== ii, :] for ii in samples]
    new_df[Symbol(string("N_", tname))] = [sum(ii[Symbol(string("N_", tname))]) for ii in sp_sam]
    new_df[Symbol(string("BF_", tname))] = new_df[Symbol(string("N_", tname))] ./ 2 ./ new_df[:N_cells]
    bs = [vcat(ii[ii[Symbol(string("locus1_int2_", tname))] .> limit, Symbol(string("locus1_int2_", tname))], ii[ii[Symbol(string("locus2_int2_", tname))] .> limit, Symbol(string("locus2_int2_", tname))]) for ii in sp_sam]
    new_df[Symbol(string("BS_mean_", tname))] = [Statistics.mean(ii) for ii in bs]
        new_df[Symbol(string("BS_median_", tname))] = [ if isempty(ii) 0 else median(ii) end for ii in bs]
    new_df[Symbol(string("BS_std_", tname))] = [Statistics.std(ii) for ii in bs]
    new_df[Symbol(string("locus_dist_mean_", tname))] = [Statistics.mean(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) for ii in sp_sam]
    new_df[Symbol(string("locus_dist_median_", tname))] = [if isempty(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) 0 else median(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) end for ii in sp_sam]
    new_df[Symbol(string("locus_dist_std_", tname))] = [std(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) for ii in sp_sam]
 
    end
    
    calculate(t, tname)

    
    new_df
    
end

function calculate_bf(t2, t3, t4)
    new_df = DataFrames.DataFrame()
    
    t2[:Sample] = t2[:Genotype] .* "_" .* string.(t2[:Timepoint])
    t3[:Sample] = t3[:Genotype] .* "_" .* string.(t3[:Timepoint])
    t4[:Sample] = t4[:Genotype] .* "_" .* string.(t4[:Timepoint])
    samples = unique(t2[:Sample])
    new_df[:Sample] = samples
    new_df[:Genotype] = [split(ii, "_")[1] for ii in samples]
    new_df[:Timepoint] = [split(ii, "_")[2] for ii in samples]
    
    function calculate(t, tname)

    new_df[:N_cells] = [sum(t[:Sample] .== ii) for ii in samples]
    sp_sam = [t[t[:Sample] .== ii, :] for ii in samples]
    new_df[Symbol(string("N_", tname))] = [sum(ii[Symbol(string("N_", tname))]) for ii in sp_sam]
    new_df[Symbol(string("BF_", tname))] = new_df[Symbol(string("N_", tname))] ./ 2 ./ new_df[:N_cells]
    bs = [vcat(ii[ii[Symbol(string("locus1_int2_", tname))] .> 0, Symbol(string("locus1_int2_", tname))], ii[ii[Symbol(string("locus2_int2_", tname))] .> 0, Symbol(string("locus2_int2_", tname))]) for ii in sp_sam]
    new_df[Symbol(string("BS_mean_", tname))] = [Statistics.mean(ii) for ii in bs]
    new_df[Symbol(string("BS_median_", tname))] = [median(ii) for ii in bs]
    new_df[Symbol(string("BS_std_", tname))] = [Statistics.std(ii) for ii in bs]
    new_df[Symbol(string("locus_dist_mean_", tname))] = [Statistics.mean(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) for ii in sp_sam]
    new_df[Symbol(string("locus_dist_median_", tname))] = [median(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) for ii in sp_sam]
    new_df[Symbol(string("locus_dist_std_", tname))] = [std(ii[ii[Symbol(string("N_", tname))] .== 2, Symbol(string("locus_dist_", tname))] ) for ii in sp_sam]
 
    end
    
    calculate(t2, :TSS2)
    calculate(t3, :TSS3)
    calculate(t4, :TSS4)
    
    new_df
    
end

function limit_tss(tss, lim)
    t = split(string(names(tss)[4]), "_")[end]
    locus1_columns = columns_containing(tss, "locus1")
    locus2_columns = columns_containing(tss, "locus2")
    
    for ii in 1:nrow(tss)
        b = Symbol("locus1_int2_"*t)
        if tss[ii, b] != "NA"
            if 0 < tss[ii, b] < lim
                 tss[ii, Symbol("N_"*t)] = tss[ii, Symbol("N_"*t)] - 1
                for col in locus1_columns
                    try
                    tss[ii, col] = 0
                    catch
                    tss[ii, col] = "0"
                    end

                end
            end
        end
        
        b = Symbol("locus2_int2_"*t)
        if tss[ii, b] != "NA"
            if 0 < tss[ii, b] < lim
                 tss[ii, Symbol("N_"*t)] = tss[ii, Symbol("N_"*t)] - 1
                for col in locus2_columns
                    try
                    tss[ii, col] = 0
                    catch
                    tss[ii, col] = "0"
                    end

                end
            end
        end
        
        
                
       
    end

    return tss
    
end

function TSS_raw_quant(t2, tss_folder, image_folder, n; xy = 0.189, zx = 0.5)
    
    images_pat = TSSs.get_image_patterns(t2)
    
    p = Progress(length(images_pat), 1)

    dfs = []

    for a in 1:length(images_pat)
        next!(p)
        

        pat = images_pat[a]

        # Get all the TSS in the images
        a = TSSs.find_outline(tss_folder, pat)
        tss = TSSs.cell_tss_dict(a)
        tss = TSSs.meassure_tss(tss, image_folder, pat, n, xy = xy, zx = zx)
    
        df = TSSs.analysis_singlechannel(tss; image = pat,  tss_name = :TSS2)
        push!(dfs, df)

    end
    
    d = TSSs.join_in_all_common_columns(dfs...)
    return d
    
end


# Quantify thresholded dots as we quantify TSSs
function quantdots(imagefolder, typefolder; r = 2)
    type_folder = readdir(typefolder, join = true)
    dots = CSV.read(type_folder[occursin.("_FISH-QUANT__threshold_spots_", type_folder)][1], DataFrames.DataFrame, skipto = 15,header=14)
    images = split_by(dots, :File)
    p = Progress(length(images))
    for imagetostudy in keys(images)
        image3D = TSSs.read_tiff_as_gray(normpath(imagesfolder,imagetostudy))
        imdots = dots[dots[!,:File].==imagetostudy, :]
        images[imagetostudy][!,:Spot_r*string(r)] = TSSs.int_brightest_pixel(image3D, imdots[:, :Y_det], imdots[:, :X_det], imdots[:, :Z_det]; radious = r)
        next!(p)
    end
    return vcat([images[ii] for ii in keys(images)]...)
end

