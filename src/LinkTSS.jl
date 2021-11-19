function split_by(df, thing)
    samples = unique(df[!,thing])
    dfs = Dict()
    for sample in samples
        f(x) = x == sample
         d = df[[f(x) for x in df[!,thing]],:]
        dfs[sample] = d
    end
    return dfs
end

function transform_in_dict_image_cell(df) 
    sp_im = split_by(df, :Image)
    for im in keys(sp_im); 
        sp_im[im] =  split_by(sp_im[im], :Cell)
    end
    
    return sp_im
end


# For 3 channels, not sure wether it works
function measure_tss_distances(t2, t3, t4; xy = 0.189, zx = 0.5)
    # Get the center of the TSS to be able to measure distances
    p2 = []; t = t2; 
    for tss in keys(t); push!(p2, [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    
    p3 = []; t = t3; 
    for tss in keys(t); push!(p3, [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    
    p4 = []; t = t4; 
    for tss in keys(t); push!(p4,  [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    ############################################################
    ## Measure all the distances between a TSS and the others ##
    ############################################################
    
    ## Initialise vatiables
    dot_distances = DataFrame(); distan = []
    dot2 = []; dot3 = []; dot4 = []
    
    i2 = 0
    for p in p2
        i2 +=1; i3 = 0
        for o in p3
            i3 +=1
            push!(dot2, collect(keys(t2))[i2])
            push!(dot3, collect(keys(t3))[i3])
            push!(dot4, 0)
            push!(distan, euclidean(p, o))
        end
        
        i4 = 0
        for q in p4
            i4 +=1
            push!(dot2, collect(keys(t2))[i2])
            push!(dot3, 0)
            push!(dot4, collect(keys(t4))[i4])
            push!(distan, euclidean(p, q))
   
        end
    end
    i3 = 0
    
    for o in p3
        i3 +=1
        i4 =0
        for q in p4
            i4 +=1
            push!(dot2, 0)
            push!(dot3, collect(keys(t3))[i3])
            push!(dot4, collect(keys(t4))[i4])
            push!(distan, euclidean(o, q))
            
        end
    end
    
    dot_distances[!,:TSS2] = dot2
    dot_distances[!,:TSS3] = dot3
    dot_distances[!,:TSS4] = dot4
    dot_distances[!,:distance] = distan
    return dot_distances
    
end

# For 2 channels
function measure_tss_distances(t2, t3; xy = 0.189, zx = 0.5)
    # Get the center of the TSS to be able to measure distances
    p2 = []; t = t2; 
    for tss in keys(t); push!(p2, [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    
    p3 = []; t = t3; 
    for tss in keys(t); push!(p3, [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    
    p4 = []; # t = t4; 
    #for tss in keys(t); push!(p4,  [t[tss][:x_], t[tss][:y_], t[tss][:z_]]); end
    ############################################################
    ## Measure all the distances between a TSS and the others ##
    ############################################################
    
    ## Initialise vatiables
    dot_distances = DataFrame(); distan = []
    dot2 = []; dot3 = []; dot4 = []
    
    i2 = 0
    for p in p2
        i2 +=1; i3 = 0
        for o in p3
            i3 +=1
            push!(dot2, collect(keys(t2))[i2])
            push!(dot3, collect(keys(t3))[i3])
            push!(dot4, 0)
            push!(distan, euclidean(p, o))
        end
        
        i4 = 0
        for q in p4
            i4 +=1
            push!(dot2, collect(keys(t2))[i2])
            push!(dot3, 0)
            push!(dot4, collect(keys(t4))[i4])
            push!(distan, euclidean(p, q))
   
        end
    end
    i3 = 0
    
    for o in p3
        i3 +=1
        i4 =0
        for q in p4
            i4 +=1
            push!(dot2, 0)
            push!(dot3, collect(keys(t3))[i3])
            push!(dot4, collect(keys(t4))[i4])
            push!(distan, euclidean(o, q))
            
        end
    end
    
    dot_distances[!,:TSS2] = dot2
    dot_distances[!,:TSS3] = dot3
    dot_distances[!,:TSS4] = dot4
    dot_distances[!,:distance] = distan
    return dot_distances
    
end




function transform_into_TSSdict(t, pat; z = true)
    tss = Dict()
    l1 = t[!,Symbol("locus1_"*pat)][1]
    l2 = t[!,Symbol("locus2_"*pat)][1]
    # If it is a real TSS
    if startswith(string(l1), "TxS")
        tss[l1] = Dict()
        tss[l1][:x_] = t[!,Symbol("locus1_x_"*pat)][1]
        tss[l1][:y_] = t[!,Symbol("locus1_y_"*pat)][1] 
        tss[l1][:z_] = t[!,Symbol("locus1_z_"*pat)][1]
        tss[l1][:int2] = t[!,Symbol("locus1_int2_"*pat)][1]
        if z == false; tss[l1][:z_] = 0. end
    
        if startswith(string(l2), "TxS")
            tss[l2] = Dict()
            tss[l2][:x_] = t[!,Symbol("locus2_x_"*pat)][1]
            tss[l2][:y_] = t[!,Symbol("locus2_y_"*pat)][1]
            tss[l2][:z_] = t[!,Symbol("locus2_z_"*pat)][1]
            tss[l2][:int2] = t[!,Symbol("locus2_int2_"*pat)][1]
            if z == false; tss[l2][:z_] = 0. end
        end
    end
    
    return tss

end

# Works for 2 channels not sure wether it works for 3
function assign_locus_from_distances(dot_distances; max_distance = 23)
    ##########################
    ## Assign TSS to locus ##
    ##########################
    
    # Initialise locus
    
    locus1 = Dict()
    locus2 = Dict()
    
    # Functions that look horrible but, actually, made my life easier
    
    function is_in_locus(n, tss)
            if n == 1; return haskey(locus1, tss) 
            elseif n == 2; return haskey(locus2, tss) 
            end
    end
        
    function is_same_tss_in_locus(n, tss, tssname)
            if n == 1; 
                return locus1[tss] ==  tssname
            elseif n == 2; 
                return locus2[tss] ==  tssname
            end
    end
    
    function add_to_locus(n, pair1, pair2)
            if n == 1; l = locus1
            elseif n ==2; l = locus2
            end
            l[pair1[1]] = pair1[2]; l[pair2[1]] = pair2[2]
            push!(locus, n)
        end
    
    function case(LOCUS, tss_type, tss_name)
        ## LEGEND ################################
        # This function returns:
        # 1 = tss_type is not in the locus
        # 2 = tss_type isin the locus and it is tss_name
        # 3 = tss_type is in the locus and it is NOT tss_name
        # 4 = tss_name  has already been assigned
        #############################################

        otherlocus = 0
        if LOCUS == 1; otherlocus = 2
        elseif LOCUS == 2; otherlocus = 1
        end
        
        if !is_in_locus(LOCUS, tss_type)
            if is_in_locus(otherlocus, tss_type)
                if is_same_tss_in_locus(otherlocus, tss_type, tss_name)
                    return 4
                else
                    return 1
                end
            else
                return 1
            end
        elseif is_in_locus(LOCUS, tss_type) && is_same_tss_in_locus(LOCUS, tss_type, tss_name)
            if is_in_locus(otherlocus, tss_type)
                if is_same_tss_in_locus(otherlocus, tss_type, tss_name)
                    return 4
                else
                    return 2
                end
            elseif !is_in_locus(otherlocus, tss_type)
                return 2
            end
        elseif is_in_locus(LOCUS, tss_type) && !is_same_tss_in_locus(LOCUS, tss_type, tss_name)
            return 3
        end
    end
    
    #dot_distances = sort!(dot_distances, cols = :distance, rev = false)
    dot_distances = sort!(dot_distances, :distance, rev = false)

    locus = []
    locus1_tss1s = []
    locus2_tss1s = []
    locus1_tss2s = []
    locus2_tss2s = []
    
        for ii in 1:nrow(dot_distances)
        tssin = []
        for col in [:TSS2, :TSS3, :TSS4]
            if dot_distances[ii, col] .!= 0
                push!(tssin, col)
            end
        end
        
        # Name variables because otherwise it will be a hell to read this code
        # Specifically the name of the tss and their type
        
        type1 = tssin[1]; type2 =  tssin[2]
        name1 =  dot_distances[ii, tssin[1]]; name2 =  dot_distances[ii, tssin[2]]
        
        # Get cases in locus 1
        locus1_tss1 = case(1, type1, name1); locus1_tss2 = case(1, type2, name2)
        push!(locus1_tss1s, locus1_tss1); push!(locus1_tss2s, locus1_tss2); 
        # Get cases in locus 2
        locus2_tss1 = case(2, type1, name1); locus2_tss2 = case(2, type2, name2)
        push!(locus2_tss1s, locus2_tss1); push!(locus2_tss2s, locus2_tss2); 
         
        if (locus1_tss1, locus1_tss2) == (1 , 1) 
            if dot_distances[ii, :distance] < max_distance
                add_to_locus(1, type1 => name1, type2 => name2)
                locus1[Symbol(string(type1, "_", type2, "_dist"))] = dot_distances[ii, :distance]
            else
                locus1[type1] = name1; locus2[type2] = name2; push!(locus, 0)
            end
            
        elseif (locus1_tss1, locus1_tss2) == (1 , 1) && dot_distances[ii, :distance] > max_distance
            locus1[type1] = name1; locus2[type2] = name2; push!(locus, 0)
            
        elseif (locus1_tss1, locus1_tss2) == (2, 1)
           
            if locus2_tss1 !=2 && locus2_tss2 !=2 
                if dot_distances[ii, :distance] < max_distance || locus2_tss1  >= 3  
                    add_to_locus(1, type1 => name1, type2 => name2)
                    locus1[Symbol(string(type1, "_", type2, "_dist"))] = dot_distances[ii, :distance]
                    
                else
                    locus1[type1] = name1; locus2[type2] = name2; push!(locus, 0)
                end
            else
                push!(locus, 0)
            end
        
                            
        
            
        elseif (locus1_tss1, locus1_tss2) == (1, 2)
            
            if locus2_tss1 !=2 && locus2_tss2 !=2
                if dot_distances[ii, :distance] < max_distance || locus2_tss2 >= 3 
                   add_to_locus(1, type1 => name1, type2 => name2)
                   locus1[Symbol(string(type1, "_", type2, "_dist"))] = dot_distances[ii, :distance]
                else
                   locus2[type1] = name1; locus1[type2] = name2; push!(locus, 0)
                end
            else
                push!(locus, 0)
            end

            
        elseif (locus1_tss1, locus1_tss2) == (2 , 2)
            add_to_locus(1, type1 => name1, type2 => name2)
            locus1[Symbol(string(type1, "_", type2, "_dist"))] = dot_distances[ii, :distance]
        
        elseif .|(locus1_tss1 == 3, locus1_tss2 == 3) &&  locus2_tss1 < 3 && locus2_tss2 < 3

            add_to_locus(2, type1 => name1, type2 => name2)
            locus2[Symbol(string(type1, "_", type2, "_dist"))] = dot_distances[ii, :distance]
        
        elseif locus1_tss1 == 3 && locus2_tss1 < 3 
            locus2[type1] = name1; push!(locus, 0)
        
        elseif locus1_tss2 == 3 && locus2_tss2 < 3 
            locus2[type2] = name2; push!(locus, 0)
            
        else
            push!(locus, 0)
        end
        
    end
    ## This helps to debug
   """ 
    dot_distances[:locus] = locus
    dot_distances[:locus1_tss1s ] =   locus1_tss1s 
    dot_distances[:locus1_tss2s ] =  locus1_tss2s
    dot_distances[:locus2_tss1s ] =  locus2_tss1s 
    dot_distances[:locus2_tss2s ] =  locus2_tss2s 
    """
    
    return locus1, locus2

end

function findpattern(rep)
    split(names(rep)[findfirst(x -> startswith(x, "locus1"), names(rep))], "_")[2]
end


# For 3 channels not sure wether it works
function link_TSSs_cell(t2, t3, t4; max_distance = 23, z = true)
    # Store number ot TSS in the cell 
    n2 =  t2[!,"N_"*findpattern(t2)][1]; n3 = t3[!,"N_"*findpattern(t3)][1]; n4 =  t4[!,"N_"*findpattern(t4)][1]
    T2 = transform_into_TSSdict(t2, findpattern(t2), z = z); T3 = transform_into_TSSdict(t3, findpattern(t3), z = z); T4 = transform_into_TSSdict(t4, findpattern(t4), z = z);
    dot_distances = measure_tss_distances(T2, T3, T4)
    locus1, locus2 = assign_locus_from_distances(dot_distances; max_distance = max_distance)
    
    # For cases when there are no other types of TSS detected
    # Note for myself: I thing this is the most complicated function I have ever wrote
    
    if n2 == 1 && n3 == 0 && n4 == 0
        locus1[:TSS2] = collect(keys(T2))[1]
    elseif n2 == 2 && n3 == 0 && n4 == 0
        locus1[:TSS2] = collect(keys(T2))[1]
        locus2[:TSS2] = collect(keys(T2))[2]
    end
    
    if n2 == 0 && n3 == 1 && n4 == 0
        locus1[:TSS3] = collect(keys(T3))[1]
    elseif n2 == 0 && n3 == 2 && n4 == 0
        locus1[:TSS3] = collect(keys(T3))[1]
        locus2[:TSS3] = collect(keys(T3))[2]
    end
    
    if n2 == 0 && n3 == 0 && n4 == 1
        locus1[:TSS4] = collect(keys(T4))[1]
    elseif n2 == 0 && n3 == 0 && n4 == 2
        locus1[:TSS4] = collect(keys(T4))[1]
        locus2[:TSS4] = collect(keys(T4))[2]
    end
    
    locus_dict = Dict(:locus1 => locus1, :locus2 => locus2)
    return locus_dict
    
end

# For 2 channels it works

function link_TSSs_cell(t2, t3; max_distance = 23, z = true)
    # Store number ot TSS in the cell 
    n2 =  t2[!,"N_"*findpattern(t2)][1]; n3 = t3[!,"N_"*findpattern(t3)][1]; n4 = 0 #n4 =  t4[!,"N_"*findpattern(t4)][1]
    T2 = transform_into_TSSdict(t2, findpattern(t2), z = z); T3 = transform_into_TSSdict(t3, findpattern(t3), z = z); #T4 = transform_into_TSSdict(t4, findpattern(t4));
     dot_distances = measure_tss_distances(T2, T3)
    nt2 = length(unique(dot_distances[dot_distances[!,:TSS2].!=0, :TSS2]))
    nt3 = length(unique(dot_distances[dot_distances[!,:TSS3].!=0, :TSS3]))
   # nt4 = length(unique(dot_distances[dot_distances[!,:TSS4].!=0, :TSS4]))
    dot_distances[!,:TSS2_n] = [if ii == 0 0 else nt2 end for ii in dot_distances[!,:TSS2]]
    dot_distances[!,:TSS3_n] = [if ii == 0 0 else nt3 end for ii in dot_distances[!,:TSS3]]
   # dot_distances[!,:TSS4_n] = [if ii == 0 0 else nt4 end for ii in dot_distances[!,:TSS4]]
    

    #return TSSs.assign_locus_from_distances(dot_distances; max_distance = max_distance)
        
    locus1, locus2 = assign_locus_from_distances(dot_distances; max_distance = max_distance)
    
    # For cases when there are no other types of TSS detected
    # Note for myself: I thing this is the most complicated function I have ever wrote
    
    if n2 == 1 && n3 == 0 && n4 == 0
        locus1[:TSS2] = collect(keys(T2))[1]
    elseif n2 == 2 && n3 == 0 && n4 == 0
        locus1[:TSS2] = collect(keys(T2))[1]
        locus2[:TSS2] = collect(keys(T2))[2]
    end
    
    if n2 == 0 && n3 == 1 && n4 == 0
        locus1[:TSS3] = collect(keys(T3))[1]
    elseif n2 == 0 && n3 == 2 && n4 == 0
        locus1[:TSS3] = collect(keys(T3))[1]
        locus2[:TSS3] = collect(keys(T3))[2]
    end
    
    if n2 == 0 && n3 == 0 && n4 == 1
        locus1[:TSS4] = collect(keys(T4))[1]
    elseif n2 == 0 && n3 == 0 && n4 == 2
        locus1[:TSS4] = collect(keys(T4))[1]
        locus2[:TSS4] = collect(keys(T4))[2]
    end
    
    locus_dict = Dict(:locus1 => locus1, :locus2 => locus2)
    return locus_dict
    
end




function frequency_tb(locus_link)
    cells = collect(keys(locus_link))
    df = DataFrame(Cell = cells)
    df[!,:locus1_tss2] = [try locus_link[cell][:locus1][:TSS2] catch; 0 end for cell in cells]
    df[!,:locus1_tss3] = [try locus_link[cell][:locus1][:TSS3] catch; 0 end for cell in cells]
    df[!,:locus1_tss4] = [try locus_link[cell][:locus1][:TSS4] catch; 0 end for cell in cells]
    
    df[!,:locus1_tss2_tss3] = [try locus_link[cell][:locus1][:TSS2_TSS3_dist] catch; "NA" end for cell in cells]
    df[!,:locus1_tss2_tss4] = [try locus_link[cell][:locus1][:TSS2_TSS4_dist] catch; "NA" end for cell in cells]
    df[!,:locus1_tss3_tss4] = [try locus_link[cell][:locus1][:TSS3_TSS4_dist] catch; "NA" end for cell in cells]
    
    df[!,:locus2_tss2] = [try locus_link[cell][:locus2][:TSS2] catch; 0 end for cell in cells]
    df[!,:locus2_tss3] = [try locus_link[cell][:locus2][:TSS3] catch; 0 end for cell in cells]
    df[!,:locus2_tss4] = [try locus_link[cell][:locus2][:TSS4] catch; 0 end for cell in cells]
    
    df[!,:locus2_tss2_tss3] = [try locus_link[cell][:locus2][:TSS2_TSS3_dist] catch; "NA" end for cell in cells]
    df[!,:locus2_tss2_tss4] = [try locus_link[cell][:locus2][:TSS2_TSS4_dist] catch; "NA" end for cell in cells]
    df[!,:locus2_tss3_tss4] = [try locus_link[cell][:locus2][:TSS3_TSS4_dist] catch; "NA" end for cell in cells]

    df[!,:TSS2_N] = [count(x -> x != 0, [df[ii,:locus1_tss2], df[ii,:locus2_tss2]]) for ii in 1:nrow(df)]
    df[!,:TSS3_N] = [count(x -> x != 0, [df[ii,:locus1_tss3], df[ii,:locus2_tss3]]) for ii in 1:nrow(df)]
    df[!,:TSS4_N] = [count(x -> x != 0, [df[ii,:locus1_tss4], df[ii,:locus2_tss4]]) for ii in 1:nrow(df)]
    df
end

function quant(tss, tssname,  cell, key)
    try
        tss[cell][tssname][key]
        catch; 0
    end
end


function add_tss_quantification(t2,t3,t4,freq_tb; quantmethod = :int2)
    cells = freq_tb[!,:Cell]
    
    freq_tb[!,:locus1_TSS2_size] = [quant(t2, freq_tb[ii, :locus1_tss2], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS2_size] = [quant(t2, freq_tb[ii, :locus2_tss2], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS3_size] = [quant(t3, freq_tb[ii, :locus1_tss3], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS3_size] = [quant(t3, freq_tb[ii, :locus2_tss3], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS4_size] = [quant(t4, freq_tb[ii, :locus1_tss4], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS4_size] = [quant(t4, freq_tb[ii, :locus2_tss4], freq_tb[ii, :Cell], quantmethod) for ii in 1:nrow(freq_tb)]

    freq_tb
       
end

function quant2(tss, tssname, cell, image,  key)
    try
        tssnamebool = tss[!,:locus1_TSS2] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus1_"*string(key)*"_TSS2")][1]
        catch;
         try 
        tssnamebool = tss[!,:locus2_TSS2] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus2_"*string(key)*"_TSS2")][1]
            catch;
            0
        end
    end
end

function quant3(tss, tssname, cell, image,  key)
    try
        tssnamebool = tss[!,:locus1_TSS3] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus1_"*string(key)*"_TSS3")][1]
        catch;
         try 
        tssnamebool = tss[!,:locus2_TSS3] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus2_"*string(key)*"_TSS3")][1]
                catch;
            0
        end
    end
end

function quant4(tss, tssname, cell, image,  key)
    try
        tssnamebool = tss[!,:locus1_TSS4] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus1_"*string(key)*"_TSS4")][1]
        catch;
        try 
        tssnamebool = tss[!,:locus2_TSS4] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("locus2_"*string(key)*"_TSS4")][1]
            catch;
            0
        end
    end
end


# For 3 channels I am not sure it works


function link_TSSs(t2, t3, t4; max_distance = 23, z = true)
    tt2 = TSSs.transform_in_dict_image_cell(t2); tt3 = TSSs.transform_in_dict_image_cell(t3); tt4 = TSSs.transform_in_dict_image_cell(t4) 
    linked_tss = Dict()
    p = Progress(nrow(t2))
    
    for im in keys(tt2)
        linked_tss[im] = Dict()
        # In case the image has been eliminated from one analysis skip it
        if haskey(tt3, im) && haskey(tt2, im)
            for cell in keys(tt2[im])
                next!(p)
                
                cell_tss2 = tt2[im][cell]; cell_tss3 = tt3[im][cell]; cell_tss4 = tt4[im][cell]
                tss2s = collect(keys(TSSs.transform_into_TSSdict(cell_tss2, findpattern(t2), z = z)))
                tss3s = collect(keys(TSSs.transform_into_TSSdict(cell_tss3, findpattern(t3), z = z)))
                tss4s = collect(keys(TSSs.transform_into_TSSdict(cell_tss4, findpattern(t4), z = z)))
                
                loc = TSSs.link_TSSs_cell(cell_tss2, cell_tss3, cell_tss4; max_distance = max_distance)
                # This bit of code was used on 20thAugust to debug the code
                """
                tss2s_assigned = [] 
                if haskey(loc[:locus1], :TSS2); push!(tss2s_assigned, loc[:locus1][:TSS2]) end
                if haskey(loc[:locus2], :TSS2); push!(tss2s_assigned, loc[:locus2][:TSS2]) end
                
                tss3s_assigned = [] 
                if haskey(loc[:locus1], :TSS3); push!(tss3s_assigned, loc[:locus1][:TSS3]) end
                if haskey(loc[:locus2], :TSS3); push!(tss3s_assigned, loc[:locus2][:TSS3]) end
                
                tss4s_assigned = [] 
                if haskey(loc[:locus1], :TSS4); push!(tss4s_assigned, loc[:locus1][:TSS4]) end
                if haskey(loc[:locus2], :TSS4); push!(tss4s_assigned, loc[:locus2][:TSS4]) end
                
                is_tss2_fine = length(intersect(tss2s_assigned, tss2s)) == length(tss2s)
                is_tss3_fine = length(intersect(tss3s_assigned, tss3s)) == length(tss3s)
                is_tss4_fine = length(intersect(tss4s_assigned, tss4s)) == length(tss4s)
                
                if is_tss2_fine && is_tss3_fine && is_tss4_fine
                    linked_tss[im][cell] = Dict()
                    linked_tss[im][cell] =  loc  
                else
                   return  tss3s_assigned, tss3s
                    #cell_tss3
                    #cell_tss2 , cell_tss3, cell_tss4
                    return loc
                
                    T2 = transform_into_TSSdict(cell_tss2, "TSS2"); T3 = transform_into_TSSdict(cell_tss3, "TSS3"); T4 = transform_into_TSSdict(cell_tss4, "TSS4");
                    return dot_distances = measure_tss_distances(T2, T3, T4)
                end
                """
                
                  linked_tss[im][cell] = Dict()
                linked_tss[im][cell] =  loc  
                

                   
            end
        end 
    end
    
    dfs = []
    
    for im in keys(linked_tss)
        fq_tb = TSSs.frequency_tb(linked_tss[im])
        fq_tb[!,:Image] = [im for ii in 1:nrow(fq_tb)]
        
        push!(dfs, fq_tb)
    end
    return fq_tb = join_in_all_common_columns(dfs...)
    
    return fq_tb = add_tss_quantification2(t2, t3, t4, fq_tb; quantmethod = :int2)
    

end


# For 2 channels it works


function link_TSSs(t2, t3; max_distance = 23, z = true)
     tt2 = TSSs.transform_in_dict_image_cell(t2); tt3 = TSSs.transform_in_dict_image_cell(t3); #tt4 = TSSs.transform_in_dict_image_cell(t4) 
    linked_tss = Dict()
    p = Progress(nrow(t2))
    
    for im in keys(tt2)
        linked_tss[im] = Dict()
        # In case the image has been eliminated from one analysis skip it
        if haskey(tt3, im) && haskey(tt2, im)
            for cell in keys(tt2[im])
                next!(p)
                
                #cell = "Cell_CP_17"
                
                cell_tss2 = tt2[im][cell]; cell_tss3 = tt3[im][cell];# cell_tss4 = tt4[im][cell]
                tss2s = collect(keys(TSSs.transform_into_TSSdict(cell_tss2, findpattern(t2), z = z)))
                tss3s = collect(keys(TSSs.transform_into_TSSdict(cell_tss3, findpattern(t3), z = z)))
                #tss4s = collect(keys(TSSs.transform_into_TSSdict(cell_tss4, findpattern(t4), z = z)))
                
                loc = TSSs.link_TSSs_cell(cell_tss2, cell_tss3; max_distance = max_distance)
             
                # This bit of code was used on 20thAugust to debug the code
                """
                tss2s_assigned = [] 
                if haskey(loc[:locus1], :TSS2); push!(tss2s_assigned, loc[:locus1][:TSS2]) end
                if haskey(loc[:locus2], :TSS2); push!(tss2s_assigned, loc[:locus2][:TSS2]) end
                
                tss3s_assigned = [] 
                if haskey(loc[:locus1], :TSS3); push!(tss3s_assigned, loc[:locus1][:TSS3]) end
                if haskey(loc[:locus2], :TSS3); push!(tss3s_assigned, loc[:locus2][:TSS3]) end
                
                tss4s_assigned = [] 
                if haskey(loc[:locus1], :TSS4); push!(tss4s_assigned, loc[:locus1][:TSS4]) end
                if haskey(loc[:locus2], :TSS4); push!(tss4s_assigned, loc[:locus2][:TSS4]) end
                
                is_tss2_fine = length(intersect(tss2s_assigned, tss2s)) == length(tss2s)
                is_tss3_fine = length(intersect(tss3s_assigned, tss3s)) == length(tss3s)
                is_tss4_fine = length(intersect(tss4s_assigned, tss4s)) == length(tss4s)
                
                if is_tss2_fine && is_tss3_fine && is_tss4_fine
                    linked_tss[im][cell] = Dict()
                    linked_tss[im][cell] =  loc  
                else
                   return  tss3s_assigned, tss3s
                    #cell_tss3
                    #cell_tss2 , cell_tss3, cell_tss4
                    return loc
                
                    T2 = transform_into_TSSdict(cell_tss2, "TSS2"); T3 = transform_into_TSSdict(cell_tss3, "TSS3"); T4 = transform_into_TSSdict(cell_tss4, "TSS4");
                    return dot_distances = measure_tss_distances(T2, T3, T4)
                end
                """
                
                  linked_tss[im][cell] = Dict()
                linked_tss[im][cell] =  loc  
                

                   
            end
        end 
    end
    
    dfs = []
    
    for im in keys(linked_tss)
        fq_tb = TSSs.frequency_tb(linked_tss[im])
        fq_tb[!,:Image] = [im for ii in 1:nrow(fq_tb)]
        
        push!(dfs, fq_tb)
    end
    return fq_tb = vcat(dfs...)    

end



function add_tss_quantification2(t2,t3,t4,freq_tb; quantmethod = :int2)
    cells = freq_tb[!,:Cell]
    
    freq_tb[!,:locus1_TSS2_size] = [quant2(t2, freq_tb[ii, :locus1_tss2], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS2_size] = [quant2(t2, freq_tb[ii, :locus2_tss2], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS3_size] = [quant3(t3, freq_tb[ii, :locus1_tss3], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS3_size] = [quant3(t3, freq_tb[ii, :locus2_tss3], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS4_size] = [quant4(t4, freq_tb[ii, :locus1_tss4], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS4_size] = [quant4(t4, freq_tb[ii, :locus2_tss4], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]

    freq_tb
       
end


function quant_r2(tss, tssname, cell, image,  key)
    try
        tssnamebool = tss[!,"locus1_"*findpattern(tss)] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("TSS1_"*string(key))][1]
        catch;
         try 
        tssnamebool = tss[!,"locus2_"*findpattern(tss)] .== tssname
        cellbool = tss[!,:Cell] .== cell
        imagebool = tss[!,:Image] .== image
        bool = tssnamebool .* cellbool .* imagebool
        tss[bool, Symbol("TSS2_"*string(key))][1]
            catch;
            0
        end
    end
end


function add_tss_quantification3(t2,t3,t4,freq_tb; quantmethod = :r2)
    cells = freq_tb[!,:Cell]
    ii = 1
    
    freq_tb[!,:locus1_TSS2_size] = [quant_r2(t2, freq_tb[ii, "locus1_tss2"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS2_size] = [quant_r2(t2, freq_tb[ii, "locus2_tss2"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS3_size] = [quant_r2(t3, freq_tb[ii, "locus1_tss3"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS3_size] = [quant_r2(t3, freq_tb[ii, "locus2_tss3"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS4_size] = [quant_r2(t4, freq_tb[ii, "locus1_tss4"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS4_size] = [quant_r2(t4, freq_tb[ii, "locus2_tss4"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]

    freq_tb
       
end



function add_tss_quantification3(t2,t3,freq_tb; quantmethod = :r2)
    cells = freq_tb[!,:Cell]
    ii = 1
    
    freq_tb[!,:locus1_TSS2_size] = [quant_r2(t2, freq_tb[ii, "locus1_tss2"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS2_size] = [quant_r2(t2, freq_tb[ii, "locus2_tss2"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    
    freq_tb[!,:locus1_TSS3_size] = [quant_r2(t3, freq_tb[ii, "locus1_tss3"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    freq_tb[!,:locus2_TSS3_size] = [quant_r2(t3, freq_tb[ii, "locus2_tss3"], freq_tb[ii, :Cell],freq_tb[ii, :Image], quantmethod) for ii in 1:nrow(freq_tb)]
    

    freq_tb
       
end

function link_TSS_GENE_ENHANCER(nascent,ehn, genefolder; suffix = "", max_distance = 5, z = false)
    tss_gene = DataFrame(CSV.File(normpath(genefolder, nascent*".csv")))
    tss_ehn = DataFrame(CSV.File(normpath(genefolder, ehn*".csv")))
    tsslinked = TSSs.link_TSSs(tss_gene, tss_ehn; max_distance = max_distance, z = z);
    tsswithquant = TSSs.add_tss_quantification3(tss_gene, tss_ehn, tsslinked; quantmethod = :r2)
    tb = fixnames_gene_enh_pairs(tsswithquant, nascent, ehn)
tb[!,:Image_Cell] = tb[!,:Image].*"__".*tb[!,:Cell]
metadata = tss_gene[!,[
        "Image_Cell", 
        "AREA_cell",
        "AREA_nuc",
         "Genotype",
        "Timepoint",
        "Rep",
        "Sample"
        ]]
if in("linkedlocus", readdir(normpath(genefolder,".."))) == false
mkdir(normpath(genefolder,"..", "linkedlocus"))
end
tb = innerjoin(metadata, tb, on = :Image_Cell)
CSV.write(normpath(genefolder, "../linkedlocus/", nascent*"__"*ehn*"__linked"*string(max_distance)*"um"*suffix*".csv"), tb)
end
