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

function x = ibvht(r,y,varargin)
% function x = ibvht(r,y,varargin)
% Inverse Balanced Vandermonde Hermitian Transform
%
% x = inv(vandermonde(r)')*y, where vandermonde(r) returns a vandermonde
% matrix with exponential series of r(k) on its rows.
% 
% Note that often, you would like to use this in combination with
% leja_indices with
%
%     N = length(x);
%     ix = leja_indices(N);
%     y(ix,1) = ibvht(r(ix),x);
%
% The above is analytically equivalent with y = ibvht(r,x), but numerically
% much more stable. Since the indices "ix" in the above example is a
% constant vector, they are calculated outside the current function.
%
% Currently two different implementations are available. The default is
% an implementation of the Björck-Pereyra algorithm from
%
% Björck, Åke, and Pereyra, Victor. "Solution of Vandermonde systems 
% of equations." Mathematics of Computation 24.112 (1970): 893-903.
%
% This implementation is invoked by default, but can also be explicitly
% accessed using the optional argument
%
%   y = ibvht(...,'method','bjorck');
%
% An alternative algorithm uses the explicit formula for the inverse from
%
% Macon, N., and A. Spitzbart. "Inverses of Vandermonde matrices." The 
% American Mathematical Monthly 65.2 (1958): 95-100.
%
% This alternative implementation can be accessed by
%
%   y = ibvht(...,'method','macon');
%
% When using the Macon implementation, an extra parameter 'threshold' can
% also be given, which indicates the regularization parameter for
% divisions (default = 2^-12). A larger threshold gives more stable but
% less accurate solutions. This parameter can be given by
%
%   y = ibvht(...,'method','Macon','threshold',2^10);
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
        x = ibvht_bjorck(r(:),y(:));
    case 1
        x = ibvht_macon(r,y,thresh);
end


end

function x = ibvht_bjorck(r,y)
N = length(r)-1;

c = zeros(N+1,2);
d = zeros(N+1,1);
c(:,1) = y;
r = conj(r); % this algorithm is for the transpose, so conjugate everything to get the hermitian

for k=0:(N-1)
    ix1 = N:-1:k+1;
    ix2 = k:-1:0;
    c(1+ix1,2) = (c(1+ix1,1) - r(1+k)*c(1+ix1-1,1));
    c(1+ix2,2) = c(1+ix2,1);
    c(:,1) = c(:,2);
end
for k=(N-1):-1:0
    ix1 = 0:k;
    ix2 = (k+1):N;
    d(1+ix1,1) = c(1+ix1,2);
    d(1+ix2,1) = c(1+ix2,2)./(r(1+ix2) - r(1+ix2-k-1));
    ix1 = [0:(k-1) N];
    ix2 = (k):(N-1);
    c(1+ix1,1) = d(1+ix1,1);
    c(1+ix2,1) = d(1+ix2,1) - d(1+ix2+1,1);
    c(:,2) = c(:,1);
end
x = c(:,1);

end



function x = ibvht_macon(r,y,thresh)
N = length(r);

% evaluate d(k) = polyeval(poly(r(h neq k)), r(k))
d = zeros(N,1);
for k=1:N
    d(k,1) = prod(1 - r(k)*r([1:(k-1) (k+1):N])');
end

% regularize d
maxd = max(abs(d));
ix = find(abs(d) < maxd*thresh);
d(ix) = maxd*thresh*(d(ix)./abs(d(ix)));

% generate A'
Ah = vandinvpolyh(r);

% multiply with A'
x = Ah*y;

% multiply y with inverse of D
x = x./conj(d);

end

function A = vandinvpolyh(r,p)
    if nargin < 2
        p = 1;
    end
    N = length(r);
    if N > 1
        N2 = ceil(N/2);
        p1 = conv(p,poly(r((N2+1):N)'));
        p2 = conv(p,poly(r((1:N2))'));
        A = [vandinvpolyh(r(1:N2),p1); vandinvpolyh(r((N2+1):N),p2)];
    else
        A = p(:)';
    end
end
