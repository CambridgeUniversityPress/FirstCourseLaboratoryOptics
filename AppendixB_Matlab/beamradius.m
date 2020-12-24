d% -------------------------------------------------------------------
% FIT FUNCTION FOR BEAMWIDTH MEASUREMENTS
% Language: Matlab and Octave
% -------------------------------------------------------------------
% Returns the field radius of a TEM_00 mode beam at any point z 
% along the optic axis. Fit parameters are beam width. The input
% arguments w0, zw, lambda, z all need to be in the same units.
% The output arguments will be in those units.
%
% SYNTAX: [w,R,zR]=beamradius([w0,zw,lambda],z);
% 
% w0 = waist size
% zw = position of waist
% lambda = wavelength
%
% w  = spot size (field radius) at z
% R  = curvature of phasefront at z
% zR = Raleigh length.
% -------------------------------------------------------------------

function [w,R,zR]=beamradius(params,z)

w0=params(1);                           % beam (field) width at waist [meters]
zw=params(2);                           % waist position [meters]
lambda=params(3);                       % wavelength [meters]
zR=pi*w0^2/lambda;                      % Raleigh length [meters]

w=w0.*sqrt(1+((z-zw)/zR).^2);           % beam width at z [meters]

if nargout>=2
    R=z.*(1+(zR./z).^2);                % beam phasefront curvature at z
end