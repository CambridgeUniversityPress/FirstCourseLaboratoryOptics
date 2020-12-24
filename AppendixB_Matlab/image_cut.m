% -------------------------------------------------------------------------------------------
% IMAGE CUT DEMONSTRATION 
% Language: Matlab or Octave
% -------------------------------------------------------------------------------------------
% This script loads an image from a file whose full or relative path is specified. The cut
% is made vertically through the center of the image. Then the vector containing
% the cut intensities is smoothed to reduce laser speckle and displayed.
% -------------------------------------------------------------------------------------------

imagefile = 'myphoto.jpg';         	% User specified path to image file
A=imread(imagefile);               	% image is read in as an array A
A=double(A)/255;                 	% doubles needed for most purposes. Also normalize.
A=mean(A,3);                        % color array has three pages (R,G,B) so average->greyscale

col = round(size(A,2)/2);          	% the column index corresponding to the image center
cut = A(:,col);                     % this is how the cut is actually taken

softcut = movmean(cut,10);          % smooth out laser speckle using a moving mean

plot(cut); hold on;                 % plot the original data
plot(softcut,'linewidth',2);        % and the smoothed version
hold off;                       

