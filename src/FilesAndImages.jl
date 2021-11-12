
#using MATLAB

function maximumrgb(rgb)
    try
    #Works for RGB
        maximum([rgb.r, rgb.g, rgb.b])
    catch
    #Works for grayscale images
        rgb
    end
end


function read_tiff_as_gray(imfile)
    im = float.(ImageMagick.load(imfile))
    if typeof(im) == Array{ColorTypes.RGB{Float32}, 3}
        im = maximumrgb.(im)
    elseif typeof(im) == Array{ColorTypes.RGB{Float64}, 3}
        im = maximumrgb.(im)
    end
    return Gray.(im)
end



function brightest_pixel_localizarion(im)
    findmax(im)[2]
end



#function avg_dot_as_Mehdi(dotfile, radious = 2)
# mat"""
#AVG_image_info = imfinfo($dotfile);
#num_images = numel(AVG_image_info);
#AVG_image = imread($dotfile);
#for k = 1:num_images
#    AVG_image(:,:,k)=imread($dotfile, k, 'Info', AVG_image_info);
#end
    
#ind_AVG_image=find(AVG_image==max(max(max(AVG_image))));
#[i1, i2, i3] = ind2sub(size(AVG_image), ind_AVG_image(1));
#r=$radious;  
#brightest_AVG_image=AVG_image(i1-r:i1+r,i2-r:i2+r,i3-r:i3+r);
#quant = sum(sum(sum(brightest_AVG_image,1),2),3);
#$quant = quant
#    """
#@mget quant
#    end



function int_brightest_pixel(im; radious = 0)
    ind = brightest_pixel_localizarion(im)
    max_indexes = Tuple(ind)
    rad = radious
    return float(sum(
            im[max_indexes[1]-rad:max_indexes[1]+rad, 
               max_indexes[2]-rad:max_indexes[2]+rad, 
               max_indexes[3]-rad:max_indexes[3]+rad]))

     
end

function get_image_patterns(imagefolder)
    files = readdir(imagefolder)
    bool = [length(split(f, r"_C[\d].|_C[\d]_")) > 1 for f in files]
    files = files[bool]
    return unique([split(f, r"_C[\d].|_C[\d]_")[1] for f in files])
end

function find_outline(folder, pattern)
    files = readdir(folder)
    bool1 = [occursin("outline", ii) for ii in files]
    bool2 = [occursin(pattern, ii) for ii in files]
    files = files[bool1.*bool2]
    if length(files) == 1 
        return normpath(folder, files[1])
    end
end

function find_image(folder, pattern, channel)
    files = readdir(folder)
    bool1 = [!occursin("filtered_batch", ii) for ii in files]
    bool2 = [occursin(pattern, ii) for ii in files]
    bool3 = [endswith(ii, "_C"*"$channel"*".tif") || endswith(ii, "_C"*"$channel"*"_MAX.tif")  for ii in files]
    files = files[bool2.*bool1.*bool3]
     if length(files) == 1 
        return normpath(folder, files[1])
    end
end

function find_filtered_image(folder, pattern, channel)
     files = readdir(folder)
    #bool1 = [!occursin("filtered_batch", ii) for ii in files]
    bool2 = [occursin(pattern, ii) for ii in files]
    bool3 = [endswith(ii, "_C"*"$channel"*"_filtered_batch.tif") for ii in files]
    files = files[bool2.*bool3]
     if length(files) == 1 
        return normpath(folder, files[1])
    end
end



 function max_projection(imagetomodify; multiply = 1)
    
    dim_im = size(imagetomodify)
    new_im = zeros(Float64, dim_im[1], dim_im[2])
    
    for ii in 1:dim_im[1]
        for jj in 1:dim_im[2]
        new_im[ii, jj] = multiply*maximum(imagetomodify[ii, jj, :])
    end
    end
    new_im
    
end
