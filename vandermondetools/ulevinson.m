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

function p = ulevinson(r)
% function p = ulevinson(r)
% Levinson-Durbin algorithm for finding the solution to
%
%  toeplitz(r) * p = ones(length(r),1)

% Bäckström, Tom, International Audio Laboratories Erlangen, 2014 (c)


r = r(:);
p = zeros(size(r));
pp = p;
pnew = p;

pp(1) = 1/r(1);
p(1:2) = [1;1]/(r(1)+r(2));

for k=3:length(r)
    q = [p(1:k-1)'*[r(2:k) r(1:k-1)] pp(1:k-2)'*[r(2:k-1) r(1:k-2)]];
    g1 = (q(4) - q(3))/(q(4)*(q(1)+q(2)) - 2*q(2)*q(3));
    g2 = (1-g1*2*q(2))/q(4);
    pnew(1:k) = g1*([p(1:k-1)' 0] + [0 p(1:k-1)']) + g2*[0 pp(1:k-2)' 0];
    pp(1:k-1) = p(1:k-1);
    p(1:k) = pnew(1:k);    
end

