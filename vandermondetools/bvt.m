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

function y = bvt(r,x)
% function y = bvt(r,x)
% Balanced Vandermonde Transform
%
% Implicit evaluation of y = vandermonde(r)*x, where vandermonde(r) is a
% function which returns exponential series of r on its rows.

N = length(r(:));
y = zeros(size(x));
for k=1:N
    xtmp = x;
    s = r(k);
    while length(xtmp) > 1
        if mod(length(xtmp),2) == 0
            xtmp = xtmp(1:2:end) + s*xtmp(2:2:end);
        else
            xtmp = [xtmp(1:2:end-1) + s*xtmp(2:2:end); xtmp(end)];
        end
        s = s*s;
    end
    y(k) = xtmp;
end
