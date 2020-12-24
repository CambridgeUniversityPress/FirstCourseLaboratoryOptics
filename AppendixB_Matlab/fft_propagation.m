% ------------------------------------------------------------------------------
% FOURIER OPTICS DEMONSTRATION
% Language: Matlab or Octave
% ------------------------------------------------------------------------------ 
% Calculates the optical intensity on a screen due to aperture diffraction. 
% The complex scalar field u(x,y) in the  source plane is a concave sperical 
% phasefront (ROC=0.5 m) passing through a circular aperture. The resulting 
%complex scalar field amplitude u'(x',y') in the "field plane" is calculated 
% via Fourier transform.
% ------------------------------------------------------------------------------

% --------------------
% Physical Parameters
% --------------------
c = 3e8;                                            % speed of light in m/s
epsilon0 = 8.854e-12;                               % vacuum permittivity in F/m
lambda = 633e-9;                                    % optical wavelength in m

% -------------
% Source plane
% -------------
xpmax=0.002;                                        % Src plane area: 4*xmax*ymax m^2
ymax=xpmax;                                         
Nx = 2^nextpow2(512);                               % #pts in source plane grid = Nx*Ny
Ny = 2^nextpow2(512);                               % (nextpow2... gives next power of two for speed)
dxp = 2*xpmax/(Nx-1);  dyp=2*ymax/(Ny-1);           % interpixel dist. in the src plane (m)
xp  = repmat( ((0:Nx-1)-floor(Nx/2))  *dxp, Ny,1);  % x' values at which to calc. source field
yp  = repmat( ((0:Ny-1)-floor(Ny/2)).'*dyp, 1,Nx);  % y' values at which to calc. source field

% -----------------------
% ABCD Matrix Components
% -----------------------
% Optical system consists of [ FREE SPACE : LENS : FREE SPACE ].  
L1 = 0.035 ;%0.1e-3;                              	% aperture-lens dist. in m
L2 = 0.05;                                          % lens-screen dist. in m
f  = -0.03;                                         % f=Inf corresponds to no lens 
M = [[1 L2];[0 1]] * [[1 0];[-1/f 1]] * [[1 L1];[0 1]]; % ABCD matrix of the system
AA = M(1,1); BB = M(1,2); CC = M(2,1); DD = M(2,2);     % The components A, B, C, D

% ---------
% Aperture
% ---------
% Field amplitude is non-zero at these values of x, y (i.e. where it passes 
% through the aperture). The apertures are defined as logical matrixes that 
% are used to index the source field distribution, i.e. Usource(~aperture)=0;
% UIsource(aperture)= <something nonzero>.

% a = 50*1e-6;                                     	% circular obstrution diam. (m)
% b = 600e-6;
% aperture = (xp+0.75*b).^2+(yp-0.35*b).^2 > (a/2)^2; % circular obstruction logical mask

a = 3000*1e-6;                                      % equil.triang. aperture side (m)
aperture = ((yp<sqrt(3)*xp+a/2/sqrt(3)) &...
            (yp<-sqrt(3)*xp+a/2/sqrt(3)) &...
            (yp>-a/2/sqrt(3)));                     % equil. triangular aperture

% a = 300e-6;                                     	% triangular diam. (m)
% b = 600e-6;
% aperture = ~(((yp-0.35*b)<sqrt(3)*(xp+0.75*b)+a/2/sqrt(3)) &...
%             ((yp-0.35*b)<-sqrt(3)*(xp+0.75*b)+a/2/sqrt(3)) &...
%             ((yp-0.35*b)>-a/2/sqrt(3)));           % equil. triangular obstruction

% -------------
% Source Field 
% -------------
% Here, the incident field is assumed to be a Gaussian beam of width "w" 
% and radius of curvature "roc". The beam is clipped by the apeture.
roc = 0.5;                                          % R.O.C. of phasefront at src plane (m)
w  = 750e-6;                                        % beam width of incident beam (m) 
I0 = 7617.5;                                        % max src plane intensity (W/m^2)
E0 = sqrt(2*I0/c/epsilon0);                         % m ax field ampl. in src plane (N/C)
k=2*pi/lambda;                                      % wave number
r=sqrt(xp.^2+yp.^2);                                % src plane coordss dist from center
usource = E0*exp(-r.^2/w^2).*exp(1i*k*r.^2/2/roc);  % field ampl. in src plane
usource(~aperture)=0;                               % field is zero except in the aperture
Isource = epsilon0*c/2*abs(usource).^2;             % Intensity in the source plane (W/m^2)

% ========================================================================================
% |+|+|+|+|  THE COMPUTATION OCCURS BETWEEN THIS LINE AND THE ONE LIKE IT BELOW  |+|+|+|+|
% ========================================================================================

% h, below is a scale factor to change from the physical units (meters) to new units in
% which all physical lengths are scaled by h=sqrt(B*lambda). In the new units, the Fresnel
% integral becomes a standard Fourier transform multiplied by a phase factor. We now scale
% all physical lengths to the new units before performing the fourier tranform. Due to the 
% limitations placed on variable names, x' in the text is the variable f here, y' is g, 
% X' is F, and Y' is G.

h = sqrt(BB*lambda);                                % scaling factor
dXp=dxp/h; dYp=dyp/h;                               % src interpixel dist in the new units
Xp = xp/h;                                          % src plane x-coords scaled to new units
Yp = yp/h;                                          % src plane y-coords scaled to new units

dX = 1/dXp/Nx;  dY=1/dYp/Ny;                        % corresponding spatial sampling interval 
                                                    % in field plane after 2 dim. FFT (fft2).
X=repmat(([0:Nx-1]-floor(Nx/2))  *dX,Ny,1);         % Field plane, x-domain (in scaled length)
Y=repmat(([0:Ny-1]-floor(Ny/2)).'*dY,1,Nx);         % Field plane, y-domain (in scaled length)
dx=dX*h;  dy=dY*h;                                 	% field plane sampling interval (in meters)        
x = X*h;  y = Y*h;                                  % Field plane, x and y-domains (in meters)

% Perform 2D FFT on and scale correctly
% -------------------------------------
ufield = ...
    -1i*exp(1i*pi*DD/BB/lambda*((x).^2+(y).^2))...  % Perform the 2D FFT on the field 
    .*fftshift( fft2( exp(1i*pi*AA*(Xp.^2+Yp.^2)).*usource )*dXp*dYp ); % FT2
Ifield = epsilon0*c/2*abs(ufield).^2;               % get the intensity

% ========================================================================================
% |+|+|+|+| CODE BELOW CHECKS AND DISPLAYS THE RESULTS |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
% ========================================================================================

% Check energy conservation
% -------------------------
inpow =  trapz(trapz(Isource))*dxp*dxp;				% integral of intensity in the src plane 
outpow = trapz(trapz(Ifield))*dx*dy;				% (total power) should equal field plane
disp(['Power in the source plane:   Pin  =  ',num2str(inpow*1000),' mW']);
disp(['Power in the field plane:   Pout  =  ',num2str(outpow*1000),' mW']);

% Display source plane intensity (Fig. 1)
% ---------------------------------------
figure(1);                                          % open a figure window
ax1 = surf(xp*1e3,yp*1e3,(Isource/1000));           % intens. (mW/mm^2) in src pl. (x,y in mm)
view(2);                                            % top-view orientation
xlabel('x (mm)');                                   % label the axes
ylabel('y (mm)');					
axis square                                         % show both axes on the same scale
axis tight                                          % minimize white space
set(ax1,'linestyle','none');                        % don't draw the axes
caxis([min([max(Isource(aperture)),...              % set the color axis
    min(Isource(aperture))]),max(Isource(aperture))]/1000); 
colormap(copper); 									% activate the color map
shading interp;										% looks more realistic
grid off;
cbar1=colorbar;
ylabel(cbar1,'Intensity (mW/mm^2)');
hold off;
title('Source Plane Irradiance');
set(gca,'fontsize',14);

% Display field plane intensity (Fig. 2)
% --------------------------------------
figure(2);
ax2 = surf(x*1e3,y*1e3,(Ifield/1000));              % plot the intensity in the field plane
view(2);                                            % can be rotated from this top-view
shading interp;
grid off;
xlabel('x (mm)');                               
ylabel('y (mm)');
axis square
axis tight                                          
set(ax2,'linestyle','none');
caxis(([min(min(Ifield)) (max(max(Ifield)))/1000]));
title('Diffracted Irradiance in the Field Plane');
colormap(copper); 
cbar2=colorbar;
ylabel(cbar2,'?Irradiance  (mW/mm^2)');
set(gca,'fontsize',14);
