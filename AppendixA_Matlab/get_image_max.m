% Gets the brightest pixel in a smoothed image
% n is the number of pixels corresponding
% to the linear distance over which averaging
% (smoothing) occurs. If the Image Processing
% toolbox is not present, use the 1D smoothing
% (currently commented out) instead of 
% the 2D Gaussian filter

function maxval = get_image_max(file,n)
    A = double(imread(file));
    A = mean(A,3);
    A = imgaussfilt(A,round(n/2));
    maxval = max(max(A));
