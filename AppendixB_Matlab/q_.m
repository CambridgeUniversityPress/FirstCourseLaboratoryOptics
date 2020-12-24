% ---------------------------------------------------------------
% GET THE Q-FACTOR OF A BEAM
% ---------------------------------------------------------------
% Returns the complex beam radius, q, of a Gaussian beam given
% the spot size, phasefront radius of curvature, R, 
% and wavelength, lambda.
%
% SYNTAX: q=q_(w,R <,lambda>);
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

function q=q_(w,R,varargin)

if nargin>=3, lambda=varargin{1}; else lambda=1064e-9; end


if R~=Inf
    q=pi*w.^2.*R./(pi*w.^2-1i.*R.*lambda);
else
    q=1i*pi*w.^2./lambda;
end