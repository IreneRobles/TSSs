
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



function int_brightest_pixel(im,x,y, z; radious = 0)
    maxinds = size(im)
    max_indexes = [x,y,z]
    
    x1 = if max_indexes[1]-radious >0 max_indexes[1]-radious else 1 end
    x2 = if max_indexes[1]+radious <= maxinds[1] max_indexes[1]+radious else maxinds[1] end
    
    y1 = if max_indexes[2]-radious >0 max_indexes[2]-radious else 1 end
    y2 = if max_indexes[2]+radious <= maxinds[2] max_indexes[2]+radious else maxinds[2] end
    
    z1 = if max_indexes[3]-radious >0 max_indexes[3]-radious else 1 end
    z2 = if max_indexes[3]+radious <= maxinds[3] max_indexes[3]+radious else maxinds[3] end
    return float(sum(im[x1:x2, y1:y2, z1:z2]))
end


function int_brightest_pixel(im,xs::Vector,ys::Vector,zs::Vector; radious = 0)
    return [convert(Float32, TSSs.int_brightest_pixel(im,xs[ii],ys[ii],zs[ii]; radious = radious)) for ii in 1:length(xs)]
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
