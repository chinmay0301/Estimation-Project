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

function r = leja_indices(N)
% function ix = leja_indices(N)
%
% Returns list of indices which can be used for Leja-ordering. This is an
% approximation which assumes that the indices are approximately evenly
% distributed over the unit circle.

r = 1:N;
r = [sub_leja(r(1:2:end)); sub_leja(r(2:2:end))];

end

function r = sub_leja(r)

if length(r) > 1
    r = [sub_leja(r(1:2:end)); sub_leja(r(2:2:end))];
end

end
