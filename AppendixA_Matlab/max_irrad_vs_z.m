% max_irrad_vs_z.m
% Author: A. Gretarsson
% 
% Plots the maximum irradiance (arb. units) in a beam versus propagation
% distance, z. Assumes the image filenames have purely numerical names corresp.
% to the location of the image taken (in whatever units used). An attempt is made
% to reduce laser speckle by smoothing the image. This requires the user to
% be judicious in the choice of nsmooth, the linear smoothing size in pixels. 
% 
% Requires: get_image_max.m


image_folder = 'sample_images';
image_extension = '.tif';
nsmooth = 32;

images = dir(image_folder);
posvals = [];
maxvals = [];

for s = 1:length(images)
    [fdir,fname,fext] = fileparts(images(s).name);
    fullpath = fullfile( ...
               images(s).folder,images(s).name);
    if strcmp(fext,image_extension)
        posvals = [posvals, ...
               str2num(fname)];
        maxvals = [maxvals, ...
               get_image_max(fullpath,nsmooth)];
    end
end

plot(posvals,maxvals,'s','linewidth',2);
grid('on');
xlabel('Position  ( mm )');
ylabel('Max. Irradiance  ( arb. units )');
shg;