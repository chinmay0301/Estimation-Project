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

function V = vandermonde(x,M)
% function V = vandermonde(x)
% Create Vandermonde matrix with elements v, that is, a matrix which has
% exponential series in its columns.
% 
% Alternative format
%   function V = vandermonde(x,M)
% Create Vandermonde matrix of size NxM (default M=N=length(v) )


N = length(x);

if nargin < 2
    M = N;
end
V = zeros(N,M);

for k=1:N
    V(k,:) = x(k).^(0:M-1);
end
