function [v, runtime] = find_vand_recursive_alg(R,N,noise_factor)
% function [v, runtime] = find_vand_recursive_alg(R,N,noise_factor)
%
% To be used as a sub-function of 'find_vand()'.
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.

R(1) = noise_factor*R(1);
hvec = R([N:-1:2 1:N],1)';

tic;
[alph0,beta0,F0]=AlgII([hvec,0]);            %% get LU factors for H
[H,h,a]=buildH(hvec,F0);                     %% get a and h_{2n}
[alpha,beta2,F,T]=AlgII([hvec,h(end)]);      %% create tridiagonal T
try
    [e2,tt,iters]=do_lr2([[0;ones(length(beta2),1)],alpha,[beta2;0]]); %get eigs
    if max(abs(1-abs(e2(:)))) > .5
        e2 = exp(i*pi*2*(0:N-1)'/N);
    end
catch
    e2 = exp(i*pi*2*(0:N-1)'/N);
end
v=e2./abs(e2);      %% =====  try artificially projecting to unit circle.
%d=FastVandermonde(v,hvec(1:N));             %% solve for d=diag(D)

runtime = toc;



function [eigs,T0,iters]=do_lr2(T0,tol,iter_max);
% function [eigs,T0,iters]=do_lr2(T0,tol,iter_max);
% Compute eigenvalues for an arbitrary tridiagonal matrix.
% Algorithm used: implicit LR Algorithm with single shifts.
%     A random shift is used if a small pivot occurs plus every 10 iterations.
% tol defaults to [-1e-6,1e-06].  If a single number, it's used for both.
%      If tol(1)<0, then use |tol| and turn off messages.
%      tol(2) is used to test convergence , tol(1) to test for small pivots.
%      Making them too different appears to be a waste of effort.
% T0 can be a square matrix, or just the 3 diagonal strips.  Example
%    A=[1 2 0 0 0;3 4 5 0 0; 0 6 7 8 0; 0 0 9 10 11; 0 0 0 12 13]
% is equivalent to supplying the matrix by diagonal strips:
%    T= [0 3 6 9 12; 1 4 7 10 13; 2 5 8 11 0]'
% iter_max = iteration limit (default 30*n)
% note: The strips are stored vertically.
%       Internally, the subsub-diagonal is also stored (used in algorithm).
% If there are 3 or 4 rows or 3 or 4 columns, and the matrix is not square
%       then assume the input matrix T0 is stored by strips.
% OUTPUTS: eigs = computed eigenvalues (unless convergence failed)
%          T0 = input matrix converted to tridiagonal strips
%          iters = number of iterations actually performed.
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.


verbose=1;
if nargin<2, tol=[];end; if isempty(tol),tol=[-1e-6,1e-06];end;
if tol(1)<0, verbose=0;tol=abs(tol);end;
if length(tol)==1,tol(2)=tol(1);end;
[n,m]=size(T0);
eigs=zeros(n,1);
retry=0;
done=0;
iters=0;

if n==m,
   T=zeros(n,4);
   T(2:end,2)=diag(T0,-1);
   T(:,3)=diag(T0,0);
   T(1:end-1,4)=diag(T0,1);
elseif m==3,
   T=zeros(n,4);
   T(:,2:4)=T0;
elseif n==3,
   T=zeros(n,4);
   T(:,2:4)=T0.';
elseif m==4;
   T=T0;
   n=size(T,1);
elseif n==4;
   T=T0.';
   n=size(T,1);
end;
if nargin<3, iter_max=[];end ; if isempty(iter_max),iter_max=30*n;end

if nargout>1, T0=T; end;

while ~done,
   iters=iters+1;
   if iters>iter_max, disp(['do_lr2 did not converge']); eigs=T; break; end;
   if n==2, 
            eig2=eig([T(n-1,3:4);T(n,2:3)]); %was: eig2=eig(T(n-1:n,n-1:n));
	    eigs(1:2)=eig2;done=1;
	    if verbose, disp([iters,n]);disp(eigs(1:2));end; break;
   elseif abs(T(n,2))<tol(2), %was: abs(T(n,n-1))<tol,
      eigs(n)=T(n,3); %was: T(n,n);
      if verbose; disp([iters,n]); disp(eigs(n)); end;
      n=n-1;
   end;
   if retry | rem(iters,10)==0,
      shift=rand(1)-.5;
      retry = 0;
   else
      eig2=eig([T(n-1,3:4);T(n,2:3)]); %was: eig(T(n-1:n,n-1:n));
      if abs(eig2(1)-T(n,3))< abs(eig2(2)-T(n,3)), shift=eig2(1);
      %was: if abs(eig2(1)-T(n,n))< abs(eig2(2)-T(n,n)), shift=eig2(1);
      else shift=eig2(2);
      end;
   end;

   TS=T;

   mult=(T(1,3)-shift)\T(2,2); %was: mult=(T(1,1)-shift)\T(2,1);
   if abs(mult)*tol(1)>1,
      retry=1;
      T=TS;
   else
       T(2,2:3)=T(2,2:3)-mult*T(1,3:4); %was: T(2,1:3)=T(2,1:3)-mult*T(1,1:3);
       %%was: but T(1,3)==0 because the supsup diag is never touched;
       T(1,3)=T(1,3)+mult*T(1,4);
       T(2,2)=T(2,2)+mult*T(2,3);
       T(3,1)=T(3,1)+mult*T(3,2);
       %was: T(1:3,1)=T(1:3,1)+mult*T(1:3,2);
    
       for i=1:n-2,
          mult=T(i+1,2)\T(i+2,1); %was: mult=T(i+1,i)\T(i+2,i);
	  if abs(mult)*tol(1)>1,retry=1;T=TS;
	          if verbose, disp('break');iters,end;
	      break;
	  end;
          T(i+2,1:3)=T(i+2,1:3)-mult*T(i+1,2:4);
          %was: T(i+2,i:i+2)=T(i+2,i:i+2)-mult*T(i+1,i:i+2);
	  high=min(n,i+3);
	  if high>=i+1,
              T(i+1,3)=T(i+1,3)+mult*T(i+1,4);
	      if high >= i+2,
                 T(i+2,2)=T(i+2,2)+mult*T(i+2,3);
		 if high >= i+3,
                    T(i+3,1)=T(i+3,1)+mult*T(i+3,2);
		 end
              end
          end
          %was: T(i+1:high,i+1)=T(i+1:high,i+1)+mult*T(i+1:high,i+2);
       end;
   end;
end;



function [H,h,a,rootsa]=buildH(hvec,F);
% function [H,h,a,rootsa]=buildH(hvec);
%  build a Hankel matrix.  If len(hvec) is odd,
%  compute next entry so that all roots are unit.
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.


n2=length(hvec);
n=ceil(n2/2);
H=hankel(hvec(1:n),hvec(n:2*n-1));
a=[];rootsa=[];
if n*2 ~= n2;
%% U-F(:,1:size(U,2))   %%    U\(d.*(U'\h))-H\h
   if nargin>1,
      U=F(:,1:n); d = diag(U);
      a=backsolv(U,(d.*forwardsolv(U',ones(n,1))));
      disp(['residual H*a-1: ', num2str(norm(H*a-1,inf))])
   else
      a=H\ones(n,1);
   end;
   if nargout>3,rootsa=roots(a); end;
   hlast= (1 - H(n,2:n)*a(1:n-1))/a(n);
   hvec(2*n)=hlast;
end;
h=hvec(n+1:2*n);
h=h(:);


 function x = forwardsolv(A,b) 
%------does not check for singularity
% function x = forwardsolv(A,b) 
% Solves an lower triangular system
% by forward-substitution. 
%----------------------------------
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.

n = size(A,1); 
x = zeros(n,1); 
for i=1:n;
    x(i) =(b(i)-A(i,1:i-1)*x(1:i-1))/A(i,i);
end 


 function x = backsolv(A,b) 
%------does not check for singularity
% function x = backsolv(A,b) 
% Solves an upper triangular system
% by back-substitution. 
%----------------------------------
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.

n = size(A,1); 
x = zeros(n,1); 
for i=n:-1:1
    x(i) =(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end 



function [alpha,beta2,F,T]=AlgII(hvec);
% function [alpha,beta2,F,T]=AlgII(hvec);
% Algorithm II from LAA paper -- Tridiagonalization Algorithm
%   D. Boley, F. Luk, and D. Vandevoorde.
%   A Fast Method to Diagonalize a Hankel Matrix. 
%   Lin. Alg. &amp; Appl. 284:41-52, 1998.
% F(:,1:n) is the U in the LU factorization of Hankel matrix (w/ no pivoting).
%
% Author/Editor: Daniel Boley.  Copyright Boley 2013.  OK to re-use for any
% non-commercial purpose.

n2=length(hvec);
n=ceil(n2/2);
if n*2 ~= n2; abort; end;

F=zeros(n,2*n);

F(1,:)=hvec(:).';

alpha=zeros(n,1);
beta2=zeros(n,1);

% Algorithm II -- Tridiagonalization Algorithm
alpha(1) = F(1,2)/F(1,1) ; 
for i  = 1:n;
    if i>1,
        alpha(i) = F(i,i+1)/F(i,i) - F(i-1,i)/F(i-1,i-1);
        beta2(i) = F(i,i)/F(i-1,i-1);
    end;
    for j = i+1:2*n-i;
        F(i+1,j) = F(i,j+1) - alpha(i)* F(i,j);
	if i>1, F(i+1,j) = F(i+1,j) - beta2(i) * F(i-1,j); end;
    end;
end;

beta2=beta2(2:end);
if nargout>3,
   T=diag(alpha,0)+diag(beta2,1)+diag(ones(length(beta2),1),-1);
end;



