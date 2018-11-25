%/*************************************************************************
%
%         (C) Copyright Fraunhofer IIS (2014)
%
% This source code is protected by copyright law and international
% treaties. This source code is made available to you subject to the terms
% and conditions of the Gratuitous Limited Non-Commercial Source Code
% Evaluation License Agreement, which you have accepted to get access to
% this source code. If you have not accepted the terms and conditions
% mentioned above, then you are NOT ALLOWED to use this source code and
% any such unauthorised use may result in severe civil and criminal
% penalties, and will be prosecuted to the maximum extent possible under law.
% The terms and conditions mentioned above can be found at
% http://www.audiolabs-erlangen.de/resources/vandermonde-tools/package
%
%**************************************************************************/

function r = aberth(c,r,itr,tol)
%About
%   Function: Estimate real polynomial's roots by Aberth's method
%   Authors: Christian Fischer Pedersen and Tom Bäckström
%Input
%   c:      Vector with polynomial coefficients in descending order (real)
%   r:      Initial root vector (real or complex)
%   itr:    Max number of iterations, e.g. 100
%   tol:    Convergence tolerance, e.g. 1e-6
%Output
%   r:    Root vector (real or complex)
% Example:
%   c = randn(1,20); r = complex(cos(c),sin(c)); r=r(1:end-1); aberth(c,r,100,1e-6)

r = r(:);
if isreal(r)
    complex(r);
end

c  = c(:);
c1 = polyder(c);
deg = length(r);
w = nan(deg,1);

for i=1:itr
    q = polyval(c,r)./polyval(c1,r);    % O(N^2)
    for k=1:deg                           % O(N^2)
        s = sum(1./(r(k) - r([1:(k-1) (k+1):deg])));
        w(k) = -q(k)./(1-q(k).*s);
    end
    r = r+w;
end

