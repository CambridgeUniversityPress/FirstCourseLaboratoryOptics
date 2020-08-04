% ---------------------------------------------------------------
% GET THE Q-FACTOR OF A BEAM
% ---------------------------------------------------------------
% Returns the q-factor of a Gaussian beam given the spot size, w,
% phasefront radius of curvature, R, and wavelength, lambda.
%
% SYNTAX: qfactor=q_(w,R <,lambda>);
%
% w      = 1/e Field radius 
% R      = Radius of curvature of phasefront
% lambda = wavelength
%
% Any one of w, R and lambda may be a vectors or scalars.
% If more than one of w, R and lambda is a vector, all 
% vectors supplied must be the same size. w, R and lambda must
% all be in the same units.
% --------------------------------------------------------------

function qfactor=q_(w,R,varargin)

if nargin>=3, lambda=varargin{1}; else lambda=1064e-9; end


if R~=Inf
    qfactor=pi*w.^2.*R./(pi*w.^2-1i.*R.*lambda);
else
    qfactor=1i*pi*w.^2./lambda;
end