function findin(a, b)
    findall((in)(b), a)
end


function columns_containing(df, containing)
    list = names(df)
    bool = [occursin(containing, string(i)) for i in list]
    list[bool]
end


function join_in_all_common_columns(dataframe1, dataframe2)
    cols1 = names(dataframe1)
    cols2 = names(dataframe2)
    common_names = findin(cols1, cols2)
    columns_to_join = cols1[common_names]
    new_dataframe = DataFrame()
    for common_name in columns_to_join
        new_dataframe[!,common_name] = vcat(dataframe1[!,common_name], dataframe2[!,common_name])
    end
    return new_dataframe
end

function join_in_all_common_columns(dfs::Array{DataFrames.DataFrame,1})
    n_df = length(dfs)
    df = dfs[1]
    for i in 2:n_df
        df = NoLongerProblems.join_in_all_common_columns(df, dfs[i])
    end
    return df
end

function join_in_all_common_columns(dfs...)
    df = dfs[1]
    
    for i in 2:length(dfs)
        df = join_in_all_common_columns(df, dfs[i])
    end
    
    return df
end

function get_locus_data(exp)
    cols1 = columns_containing(exp, "locus1")
    push!(cols1, :Cell, :Sample, :Genotype, :Timepoint, :Image)
    locus1 = exp[:, cols1]
    cols2 = columns_containing(exp, "locus2")
    push!(cols2, :Cell, :Sample, :Genotype, :Timepoint, :Image)
    locus2 = exp[:, cols2]
    for col in cols2
        rename!(locus2, col => Symbol(string(replace(string(col), "locus2" => "locus1"))))
    end
    return join_in_all_common_columns(locus1, locus2)
end

function get_locus_data(exps...)
    dfs = []
    for ii in 1:length(exps)
        df = get_locus_data(exps[ii]); df[:Rep] = [ii for a in 1:nrow(df)]
        push!(dfs, df)
    end
    join_in_all_common_columns(dfs...)
end
