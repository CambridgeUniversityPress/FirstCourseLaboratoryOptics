
% -------------------------------------------------------------------------------------------
% IMAGE PROCESSING DEMONSTRATION -- REMOVE LASER SPECKLE
% Language: Matlab or Octave
% -------------------------------------------------------------------------------------------
% Performs a 2D fourier transform on the image in "imagefile", then sets to zero, all spatial 
% frequency bins with amplitude below the value given in  "threshold". The result is then 
% fourier transformed back to the spatial domain to produce the despeckled image. The beam 
% occupies a small part of the freq. domain image; speckle occupies the rest. Speckle has low 
% power per frequency bin. Suppressing all bins with low power therefore gets rid of speckle.
%
% SYNTAX:  I = imdespeckle(imagefile, threshold);
%
% NOTES: - threshold = 1 is usually a good starting point.
%        - The image is assumed to be a monochrome image. If it is not, then it is converted
%          to a monochrom image. 
%
% INPUT ARGUMENTS
% ---------------
% imagefile:    Full path and filename to the image (any format readable by "imread").
% threshold:    Frequency components with log10(mag) below threshold are discarded.
% -------------------------------------------------------------------------------------------

function despekld_image =imdespeckle(imagefile, threshold)

% The following lines read in the figure data, convert the color image to intensity and 
% convert the uint8 data type of the figure to the double (real numbers) data type used 
% for computation.
data = imread(imagefile);
if size(data,3)>1, data=rgb2gray(data); end
data = double(data);

% Perform the 2D numerical fourier transform and scale it correctly. The result is a
% picture of the image in "frequency space" (spatial frequency, that is).
N1 = size(data,1); N2=size(data,2);
F = fftshift(fft2(data)/sqrt(N1*N2));

% Threshold the fourier transformed image
pixels_below_threshold = log10(abs(F))<threshold;   % logical mask for pixels -> 0
Fthresh=F;                                          % start unthresholded
Fthresh(pixels_below_threshold)=0;                  % set pixels below threshold to 0                   

% Finally, perform the inverse transform on the thresholded data to get back
% to position space. (I.e. to get back our image.).
despekld_image = abs(ifft2(Fthresh)*sqrt(N1*N2));

% The plotting is done by the surf command. (There are also numerous figure options 
% called using commands like: view, set, axis, colorbar, and caxis. But these are
% only to make the picture prettier, easier to rotate in 3D, etc.)

figure(1);
h1=surf(abs(despekld_image));
set(h1,'linestyle','none');
colormap(bone); 
view(2);
set(gca,'color',[1 1 1]*0.3,'plotboxaspectratio',[N2/N1,1,1],...
    'ydir','reverse','xaxislocation','top','fontname','fixedwidth');
axis tight; 
axlims = axis;
title('Despeckled Image','fontsize',9,'fontname','fixedwidth');

figure(2);
h2=surf(log10(abs(Fthresh)));  set(h2,'linestyle','none');
view(2);
set(gca,'color',[1 1 1]*0.3,'plotboxaspectratio',[N2/N1,1,1],...
    'ydir','reverse','xaxislocation','top','fontname','fixedwidth');
axis square; 
axis(axlims); 
colorbar;
title('2D FFT, Thresholded & Smoothed','fontsize',9,'fontname','fixedwidth');