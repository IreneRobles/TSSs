
function cell_tss_dict(outiline_file)
    f = DelimitedFiles.readdlm(outiline_file)
    cell_starts = findall(x -> "CELL_START" == x, f[:, 1])
    cell_ends = findall(x -> "CELL_END" == x, f[:, 1])
    
    cell_dict = Dict{String,Dict{String,Dict}}()
    
    for a in 1:length(cell_starts)
        
         cell = if a != length(cell_starts) f[cell_starts[a]:cell_starts[a+1], :] else f[cell_starts[a]:end, :] end
         cellname = cell[1, 2]

         cell_tsss = findall(x -> "TxSite_START" == x, cell[:, 1])
        
        # Get TSS coordinates in each cell
        tsss = Dict()
        for a in 1:length(cell_tsss)
            tss = Dict{Symbol,Array{Float64, 1}}()
            tss_name = cell[cell_tsss[a], 2]
            #There was a typo fixed 20th April 2021
            tss[:x] = cell[cell_tsss[a]+2, 2:5]
            tss[:center_x] = [Statistics.mean(tss[:x])]
            tss[:y] = cell[cell_tsss[a]+1, 2:5]
            tss[:center_y] = [Statistics.mean(tss[:y])]
            tsss[tss_name] = tss
        end
        

        cell_dict[cellname] = tsss
    
    end
    
    cell_dict
end

function cell_outline_dict(outiline_file)
    f = DelimitedFiles.readdlm(outiline_file)
    cell_starts = findall(x -> "CELL_START" == x, f[:, 1])
    cell_ends = findall(x -> "CELL_END" == x, f[:, 1])
    
    cell_dict = Dict{String,Dict{String,Array{Int64,1}}}()
    
    for a in 1:length(cell_starts)
        
         cell = if a != length(cell_starts) f[cell_starts[a]:cell_starts[a+1], :] else f[cell_starts[a]:end, :] end
         cellname = cell[1, 2]
        
        #Outlines for cells

        pos = cell[2, 2:end]; bool = [typeof(ii) == Int64 for ii in pos]; pos = [ii for ii in pos[bool]]
        x = pos
        pos = cell[3, 2:end]; bool = [typeof(ii) == Int64 for ii in pos]; pos = [ii for ii in pos[bool]]
        y = pos
        
        #Outlines for nuclei
        
        pos = cell[7, 2:end]; bool = [typeof(ii) == Int64 for ii in pos]; pos = [ii for ii in pos[bool]]
        x_nu = pos
        pos = cell[8, 2:end]; bool = [typeof(ii) == Int64 for ii in pos]; pos = [ii for ii in pos[bool]]
        y_nu = pos
        
        cell_dict[cellname] = Dict("Cell_x" => x, "Cell_y" => y, "Nuc_x" => x_nu, "Nuc_y" => y_nu)
        
        
    
    end
    
    cell_dict
end
