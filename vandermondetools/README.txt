/*************************************************************************

         (C) Copyright Fraunhofer IIS (2014)

This source code is protected by copyright law and international
treaties. This source code is made available to you subject to the terms
and conditions of the Gratuitous Limited Non-Commercial Source Code
Evaluation License Agreement, which you have accepted to get access to
this source code. If you have not accepted the terms and conditions
mentioned above, then you are NOT ALLOWED to use this source code and
any such unauthorised use may result in severe civil and criminal
penalties, and will be prosecuted to the maximum extent possible under law.
The terms and conditions mentioned above can be found at
http://www.audiolabs-erlangen.de/resources/vandermonde-tools/package

**************************************************************************/



--------------------------------------------------------------------------
Documentation for Vandermonde tools

Author: Tom Bäckström
Includes contributions of prof. Daniel Boley in function find_vand_recursive_alg.m
which comes with a separate licence.


Abstract:
This is a collection of tools related to the Vandermonde factorization of 
Toeplitz matrices and its applications. This includes functions for finding
the factorization for a given Toeplitz matrix as well as applying a 
Vandermonde transform, the Hermitian Vandermonde transform and their
inverses in an efficient manner. 


--------------------------------------------------------------------------
Basics

The command 

    V = vandermonde(v);

creates a Vandermonde matrix such that the k'th row has an exponential
sequence of v(k), that is, V(k,:) = v(k).^(0:N-1). The variant, 
vandermonde_fast(), uses a faster approach to create the matrix which can 
lead to inaccuracies when the matrix is very large.


--------------------------------------------------------------------------
Transforms

Using V = vandermonde(r), the Vandermonde transform can then be expressed 
in four ways:

    y = V*x        is equivalent with   y = bvt(v,x)
    y = V'*x       is equivalent with   y = bvht(v,x)
    x = inv(V)*y   is equivalent with   x = ibvt(v,y)
    x = inv(V')*y  is equivalent with   x = ibvht(v,y)

However, such application of the inverses is numerically unstable for large
matrices. For improved stability, you can use Leja-ordering instead by

    N = length(x);
    ix = leja_indices(N);

whereby

    x = inv(V)*y   is equivalent with   x = ibvt(r(ix),y(ix))
    x = inv(V')*y  is equivalent with   x(ix) = ibvht(r(ix),y)

The forward transforms are numerically stable also without Leja-ordering.


--------------------------------------------------------------------------
Constructing a Vandermonde factorization

Suppose we have a signal x for which we want to find a frequency-domain 
representation where components are uncorrelated. The autocorrelation of x
then contains the correlation-information, whereby we start by calculation
of the autocorrelation

    r = xcorr(x,N);
    R = toeplitz(r(1+N+(0:N-1)));

The Vandermonde transform which decorrelates the signal can then be created
by

    v = find_vand('xcorr',R);

That this is indeed a factorization can be verified through

    V = vandermonde(v);
    D = V'*R*V;

where D should be diagonal. 

To decorrelate the signal, we can then use

   y = ibvht(v,x);

or for large matrices

   y(ix,1) = ibvht(v(ix),x);

Here, processing can be applied on y to obtain a modified signal yh and 
afterwards we can apply the opposite transform

   xh = bvht(v,yh);


--------------------------------------------------------------------------
