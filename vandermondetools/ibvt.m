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

function x = ibvt(r,y,varargin)
% function x = ibvt(r,y,varargin)
% Inverse Balanced Vandermonde Transform, that is, find the solution to
% x = inv(vandermonde(r))*y, where vandermonde(r) returns a vandermonde
% matrix with exponential series of r(k) on its rows.
%
% Currently two different implementations are available. The default is
% an implementation of the Bj?rck-Pereyra algorithm from
%
% Björck, Åke, and Pereyra, Victor. "Solution of Vandermonde systems 
% of equations." Mathematics of Computation 24.112 (1970): 893-903.
%
% This implementation is invoked by default, but can also be explicitly
% accessed using the optional argument
%
%   y = ibvt(...,'method','bjorck');
%
% An alternative algorithm uses the explicit formula for the inverse from
%
% Macon, N., and A. Spitzbart. "Inverses of Vandermonde matrices." The 
% American Mathematical Monthly 65.2 (1958): 95-100.
%
% This alternative implementation can be accessed by
%
%   y = ibvt(...,'method','macon');
%
% When using the Macon implementation, an extra parameter 'threshold' can
% also be given, which indicates the regularization parameter for
% divisions (default = 2^-12). A larger threshold gives more stable but
% less accurate solutions. This parameter can be given by
%
%   y = ibvt(...,'method','Macon','threshold',2^10);
%

method = 0;
if nargin > 2
    ix = 1;
    switch(lower(varargin{ix}))
        case 'method'
            switch(lower(varargin{ix+1}))
                case 'bjorck'
                    method = 0;
                case 'macon'
                    method = 1;
                    thresh = 2^-12;
                otherwise
                    error(['Unknown method ''' varargin{ix+1} '''.']);
            end
            ix = ix+2;
        case 'threshold'
            thresh = varargin{ix+1};
            ix = ix+2;
        otherwise
            error(['Unknown parameter ''' varargin{ix} '''.']);
    end
end

switch method
    case 0
        x = ibvt_bjorck(r(:),y(:));
    case 1
        x = ibvt_macon(r,y,thresh);
end

end




function x = ibvt_bjorck(r,y)

N = length(r);

c = zeros(N,2); % this probably would be possible also with a Nx1 buffer, but it's easier this way
c(:,1) = y(:);
for k=1:(N-1)
    ix1 = k:-1:1;
    ix2 = N:-1:(k+1);
    c(ix1,2) = c(ix1,1); 
    c(ix2,2) = (c(ix2,1) - c(ix2-1,1)) ./ (r(ix2) - r(ix2-k));
    c(:,1) = c(:,2);
end
for k=(N-1):-1:1
    ix1 = [1:(k-1) N];
    ix2 =  k:(N-1);
    c(ix1,1) = c(ix1,2);
    c(ix2,1) = c(ix2,2) - r(k)*c(ix2+1,2);
    c(:,2) = c(:,1);
end
x = c(:,1);


end



function x = ibvt_macon(r,y,thresh)
N = length(r);

% evaluate d(k) = polyeval(poly(r(h neq k)), r(k))
d = zeros(N,1);
for k=1:N
    %d(k) = 1;
    %for h=[1:(k-1) (k+1):N]
    %    d(k) = d(k)*(1 - r(k)*r(h)');
    %end 
    d(k,1) = prod(1 - r(k)*r([1:(k-1) (k+1):N])');
end

% regularize d
maxd = max(abs(d));
ix = find(abs(d) < maxd*thresh);
d(ix) = maxd*thresh*(d(ix)./abs(d(ix)));


% multiply y with inverse of D
yid = y./d;

% generate A
A = vandinvpoly(r);

% multiply with A
x = A*yid;


end


function A = vandinvpoly(r,p)
    if nargin < 2
        p = 1;
    end
    N = length(r);
    if N > 1
        N2 = ceil(N/2);
        p1 = conv(p,poly(r((N2+1):N)'));
        p2 = conv(p,poly(r((1:N2))'));
        A = [vandinvpoly(r(1:N2),p1) vandinvpoly(r((N2+1):N),p2)];
    else
        A = p(:);
    end
end
