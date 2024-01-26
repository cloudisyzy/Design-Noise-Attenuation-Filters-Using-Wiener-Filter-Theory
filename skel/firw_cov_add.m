% [SigmaYY,SigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noise, N)
%	
%	A		- AR model for the signal x(n), A(q)x(n)=w(n)
%	sigma2		- E[w(n)*w(n)]
%	Anoise		- AR model for the noise v(n), Anoise(q)v(n)=e(n)
%	sigma2noise	- E[e(n)*e(n)]
%	N  	- Length of Y(n), or r_yy(k_max), r is autocovariance.
%	
% 	SigmaYY		- E[Y(n) (Y(n))']
%	SigmaYx		- E[Y(n) x(n)]
%
%  firw_cov_add: Calculate covariance and cross-covariance for
%     Y(n)=[y(n), y(n-1),...,y(n-N+1)]' where y(n)=x(n)+v(n)
%     
%     Author: Ziyue Yang
%     Date: 2024.01.23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SigmaYY,SigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noise, N)

    addpath('./mfiles');

%     r_x, r_v are autocovariance of X and V, r_x=E[X(n)x(n)],
%     r_v=E[V(n)v(n)], N*1 column vector,
%     Length=N,r(0),...,r(N-1).
    r_x = ar2cov(A, sigma2, N-1);
    r_v = ar2cov(Anoise, sigma2noise, N-1);

%     E[Y(n) x(n)]=E[(X(n)+V(n))x(n)]=E[X(n)x(n)], N*1 row vector.
    SigmaYx = r_x;

%     E[Y(n) (Y(n))']=E[X(n)X(n)']+E[V(n)V(n)'], N*N matrix.
    SigmaYY = toeplitz(r_x + r_v);
    
end