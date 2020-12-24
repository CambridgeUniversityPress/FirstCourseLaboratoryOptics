% ---------------------------------------------------------------
% GET BEAM WIDTH AND ROC FROM Q
% ---------------------------------------------------------------
%?Returns the phasefront radius of curvature, R, and the beam 
% width, w, of a Gaussian beam. Accepts the complex beam radius,
% q, and the wavelength
%
% SYNTAX: [R <,w>]=R_(q <,lam>);   
%            <...> indicates optional arguments
%
% q      = q-factor of the beam at the position where R and w are to
%          be found. q can be a vector
% lam    = wavelength. Can be a vector or scalar.
%?w      = beam radius
% R      = beam phasefront curvature
%
% If both q and lam are vectors, they must be the same size.
% If w is requested as an output, lam should be supplied.
% ---------------------------------------------------------------

function [R,w]=R_(q,varargin)

if nargin>=2, lam=varargin{1}; else lam=1064e-9; end

w=sqrt(lam/pi .* imag(q).*(1+real(q).^2./imag(q).^2));
R=real(q).*(1+imag(q).^2./real(q).^2).*ones(size(w));