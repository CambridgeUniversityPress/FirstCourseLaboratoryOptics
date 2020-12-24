% --------------------------------------------------------------------
% PROPAGATE A GAUSSIAN BEAM 
% Language: Matlab and Octave
% --------------------------------------------------------------------
% Propagates a Gaussian beam (TEM_nm) with complex radius of curvature q1
% and amplitude factor p1 (optional), according to the ABCD matrix
% supplied.
%
% Returns the new complex beam radius q=(A*q1+B)/(C*q1+D) and the 
% new amplitude factor p = 1/(A+B/q1)^(1+n+m) by which the field is
% multiplied. If q1 is a vector q and p will be vectors of the same size.
%
% SYNTAX: [q,p]=prop(q1,abcd <,[n,m],p1>);
%           <...> indicates optional arguments
%
% For a Hermite Gaussian n,m are the mode designators.
% --------------------------------------------------------------------
function [q,p]=prop(q1,abcd,varargin)

if ( nargin>=3 && ~isempty(varargin{1}) ), mode=varargin{1}; else mode=[0,0]; end
if nargin>=4, p1=varargin{2}; else p1=ones(size(q1)); end

A=abcd(1,1);
B=abcd(1,2);
C=abcd(2,1);
D=abcd(2,2);

n=mode(1);
m=mode(2);

q = (A*q1 + B)./(C*q1 + D);
p = p1.*exp(1i*angle(1./(A+B./q1).^(1+n+m)));