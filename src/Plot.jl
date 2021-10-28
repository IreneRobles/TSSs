


function draw_cell(cell_ouline_dict, cellname; cell_color = "black", nu_color = "gray", showcelllabel = false)
    outline = cell_ouline_dict[cellname]
    #Draw nucleous
    plot(outline["Nuc_x"], outline["Nuc_y"], c = nu_color)
    #Draw cell
    plot(outline["Cell_x"], outline["Cell_y"], c = cell_color)
    if showcelllabel
    celln = split(cellname, "_")[end]
    miny = median(outline["Cell_y"])
    minx = median(outline["Cell_x"])
    plt.annotate(celln, xy = [minx, miny], fontsize = 12)
    end
end


function draw_TSS(cell_ouline_dict, cellname; tss_color = "purple")
    outline = cell_ouline_dict[cellname]
    if length(keys(outline)) > 0
        for tss in keys(outline)
             plot(
                vcat(outline[tss][:y], outline[tss][:y][1]),
                vcat(outline[tss][:x], outline[tss][:x][1]),
                c = tss_color
            )
        end
    end

end

