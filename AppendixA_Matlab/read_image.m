A = imread('myphoto.jpg');
A = double(A);
A=mean(A,3);
pcolor(A);
shading('flat'); colormap('bone');