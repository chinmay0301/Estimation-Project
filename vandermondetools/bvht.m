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

function y = bvht(r,x)
% function y = bvht(r,x)
% Balanced Vandermonde Hermitian Transform
%
% Implicit evaluation of y = vandermonde(r)'*x, where vandermonde(r)' is a
% function which returns exponential series of r on its columns.

N = length(r(:));
r = r(:); % make it a column vector
x = conj(x(:)); % make it a column vector

y = zeros(size(x));
xtemp = x;

y(1) = sum(xtemp)';
for k=2:N
    xtemp = xtemp.*r;
    y(k) = sum(xtemp)';
end
