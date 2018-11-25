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

function [v,d,runtime] = find_vand(varargin)
% function [v,d,runtime] = find_vand(varargin)
% Determine Vandermonde Factorization of a real values signal
%
% Example:
%
%    r = xcorr(x,N);
%    R = toeplitz(r(1+N+(0:N-1)));
%    [v,d] = find_vand('xcorr',R);
%    V = vandermonde(v);
%    D = diag(d);
%
% If you are too lazy to calculate the autocorrelation, you can give the
% input signal directly
%
%    v = find_vand('signal',x);
%

% Tom B?ckstr?m, 2013

% initialize
varix = 1;
method = 'fftescalation'; % default approach
noise_factor = 1;
oversampling = 1;
refinement = [];
stabilization = [];
iterations = 3;
while nargin >= varix
    switch lower(varargin{varix})
        case 'fftlen'   % fft length
            M = varargin{varix+1};  
            varix = varix+2;
        case 'signal'
            N = length(varargin{varix+1});
            if ~exist('M','var')
                M = N;
            end
            R = xcorr(varargin{varix+1},N);
            R = R(1+N+(0:N-1));
            R = toeplitz(R);
            varix = varix+2;
        case 'oversampling'
            oversampling = varargin{varix+1};
            varix = varix+2;
        case 'xcorr'
            R = varargin{varix+1};
            if min(size(R))==1
                R = toeplitz(R);
            end
            N = size(R,1);
            varix = varix+2;
        case 'method'
            method = varargin{varix+1};
            varix = varix+2;
        case 'noisefactor'
            noise_factor = varargin{varix+1};
            varix = varix+2;
        case 'refinement'
            refinement = varargin{varix+1};
            varix = varix+2;
        case 'stabilization'
            stabilization = varargin{varix+1};
            varix = varix+2;            
        case 'iterations'
            iterations = varargin{varix+1};
            varix = varix+2;  
        case {'levinsonmethod', 'levinson'}
            levinson_method = varargin{varix+1};
            varix = varix+2;
        case 'cond'
            c = varargin{varix+1};
            varix = varix+2;
            xi = (1+cos(pi/N))/(1 - cos(pi/N));
            M = ceil(pi/acos((xi*c - 1)/(xi*c + 1)));            
        otherwise
            error(['Unknown input ''' varargin{varix} '''.']);
    end
end
if isempty(refinement)
    if strcmp(method,'fft') || strcmp(method,'fftescalation')
        refinement = 'durand';
    else
        refinement = 'none';
    end
end
if (strcmp(method,'fft') || strcmp(method,'fftescalation')) && ~exist('M','var')
    M = 2*N;
end
if (R(1) == 0) || (size(R,1) ~= sum(isfinite(R(:,1))))
    % if input is zeros or nan give still a sane output
    method = 'uniform';
    R = eye(size(R,1));
end
if strcmp(method,'recursive_alg')
    levinson_method = 'none';
end

if ~exist('levinson_method','var')
    switch lower(method)
        case 'roots'
            levinson_method = 'inverse';
        case {'fft', 'fftescalation'}
            levinson_method = 'levinson';
        otherwise
            levinson_method = 'ulevinson';
    end
end

switch lower(levinson_method)
    case 'ulevinson'
        if sum(abs(imag(R(:)))) > 0
            error('ulevinson() implementation does not support complex-valued Toeplitz matrices.');
        end
        ar = ulevinson(R(1,:));
        a = conv(ar,[1 -1]);
    case 'levinson'
        a = levinson(R(1,:));
        a = [a'; 0] + conj(flipud([a'; 0]));
    case 'inverse'
        a = R\[1;zeros(N-1,1)];
        a = [a; 0] + conj(flipud([a; 0]));
    case 'none'
    otherwise
        error(['Unknown levinson-method: ' levinson_method])
end

switch lower(method)
    case 'roots'

        try
            tic;

            v = roots(a);
            if (max(abs(v)) > 1.1) || (min(abs(v)) < .9)
                warning('Non-unit circle zeros')
            end
            v = conj(v)./abs(v);
            ang = angle(v); ix = find(ang < 0); ang(ix) = ang(ix)+2*pi;
            [foo,ix] = sort(ang); v = v(ix);        

            runtime = toc;



        end
    case 'recursive_alg'

        [v, runtime] = find_vand_recursive_alg(R,N,noise_factor);
        
    case 'fft'
        tic;
        
        % default value for window length
        if ~exist('M','var')
            M = 10*N;       % make it quite long (less would probably be enough, but this is safer)
        end
        
        plen = length(a);
        if mod(plen,2)
            px = [a(((plen+1)/2):plen); zeros(oversampling*M-plen,1); a(1:(plen-1)/2)];
            if abs(a(1)- a(plen)) < abs(a(1)+a(plen))
                P = real(fft(px)); 
            else
                P = imag(fft(px));
            end
            P = [P; P(1)];
        else
            P = fft(a,oversampling*M).*exp(i*(0:oversampling*M-1)'*(plen-1)*pi/(oversampling*M));
            if abs(a(1)- a(plen)) < abs(a(1)+a(plen))
                P = real(P); 
            else
                P = imag(P);
            end
            P = [P;-P(1)];
        end
        ix = find(diff(sign(P)) ~= 0);  % root intervals
        if (ix(1) == 1) && (ix(end)==length(P)-1)
            ix(end) = [];
        end
        
        f = ix -1 - P(ix)./(P(ix+1) - P(ix));  % linear interpolation
        while length(f) > N
            [~,qix] = min(abs(diff(f)));
            f(qix) = [];
        end
        
        v = exp(i*2*pi*f/length(P)); % generate the complex values
        
        if length(v) < N
            warning('Incomplete decomposition');
            v = [v; eps*randn(N-length(v),1)];
        end
        
        runtime = toc;
        
    case 'fftescalation'
        
        maxmem = 100000000;
        tic;
        
        
        if strcmp(stabilization,'simple')
            error('Stabilization deactivated because it has not been tested.')
            s = stabilizing_filter(R(:,1));
            p2 = conv(s,p2);
            plen = plen+length(s)-1;
        end
        
        plen = length(a);        
        if mod(plen,2)
            px = [a(((plen+1)/2):plen); zeros(oversampling*M-plen,1); a(1:(plen-1)/2)];
            if abs(a(1)- a(plen)) < abs(a(1)+a(plen))
                partfn = @(x) (real(x));
            else
                partfn = @(x) (imag(x));
            end
            P = partfn(fft(px)); 
            Pw = P(1);
        else
            P = fft(a,oversampling*M).*exp(i*(0:oversampling*M-1)'*(plen-1)*pi/(oversampling*M));
            if abs(a(1)- a(plen)) < abs(a(1)+a(plen))
                partfn = @(x) (real(x));
            else
                partfn = @(x) (imag(x));
            end
            P = partfn(P); 
            Pw = -P(1);
        end
        ix = find(diff(sign([P;Pw])) ~= 0);  % root intervals
        
        % We want to make sure that bins are small enough such that zeros
        % are not only found, but such that they become accurate. We will
        % therefore make the FFT-length recursively longer until all zeros
        % are found, then add a couple more recursions for added accuracy.
        % This way the root-refinement algorithms converge faster.
        it = 0;
        
        while (it < 3) && (length(P) < maxmem.MaxPossibleArrayBytes/8) && (length(P) < (2^12)*N)
            if mod(plen,2)
                px = [a(((plen+1)/2):plen); zeros(M-length(a),1); a(1:((plen-1)/2))];
                Pb = fft(px .* exp(i*pi*[-(0:M/2) ((M/2)-1:-1:1)]'/M));
                Pb = partfn(Pb);
            else
                Pb = fft(a.*exp(i*pi*-(0:plen-1)'/M),oversampling*M).*exp(i*pi*(N/2+ (0:oversampling*M-1)'*N/(oversampling))/M);
                Pb = partfn(Pb);
            end
            Px = zeros(length(P)*2,1);
            Px(1:2:end) = P;
            Px(2:2:end) = Pb;
            P = Px;
            
            ix = find(diff(sign([P;Pw])) ~= 0);  % root intervals
            if length(ix) >= N % check if all roots found
                it = it+1;
            end
            M = M*2;
        end

        if ix(end) >= length(P)
            ix(end) = [];
            if ix(1) ~= 1
                ix(end+1) = 1;
            end
        end
        
        f = ix -1 - P(ix)./(P(ix+1) - P(ix));  % linear interpolation
        while length(f) > N
            [~,qix] = min(abs(diff(f)));
            f(qix) = [];
        end
        
        v = exp(i*2*pi*f/length(P)); % generate the complex values
        
        if mod(N,2)~=0
            [foo,ix] = min(abs(v+1));
            v(ix) = [];
        end
        if length(v) < N
            v = [v; eps*randn(N-length(v),1)];
        end
        
        runtime = toc;
        
        
    case 'uniform'
        runtime = 0;
        v = exp(i*2*pi*(0:N-1)/N);
        if ~strcmp(lower(refinement),'none')
            p = ulevinson(R(:,1)); p = p/p(1);     
        end
        
    otherwise
        error(['Unknown method ''' method ''' for finding Vandermonde factorization.']);
end

tic
switch lower(refinement)
    case 'none'
    case 'aberth'
        if iterations > 0
            v = aberth(conj(a),conj(v),iterations,10e-6);
        end
    case 'laguerre'
        if iterations > 0
            v = laguerre(conj(a),conj(v),iterations,10e-6);
        end
    case 'newton'
        if iterations > 0
            v = newton(conj(a),conj(v),iterations,10e-6);
        end
    case {'durand_kerner', 'durand'}
        if iterations > 0              
           v = durand_kerner(conj(a),conj(v),iterations,10e-6); 
        end
               
    otherwise
        error(['Unknown refinement method ''' refinement '''.']);
end
runtime = toc+runtime;

if nargout > 1
    if R(1) ~= 0
        d = vandermonde_fast(v)'\R(:,1);
    else
        % if input is zeros give still a sane output
        d = ones(size(v));
    end
end














