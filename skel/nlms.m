function [thetahat,xhat]=nlms(x,y,N,muu)

% [thetahat,xhat]=nlms(x,y,N,muu)
%
%	x			- Data sequence
%	y			- Data sequence
%	N			- Model order 
%	muu			- Step size
%	thetahat		- Matrix with estimates of theta. 
%				  Row n corresponds to the estimate thetahat(n)'
%	xhat			- Estimate of x
%
%
%
%  nlms: The Normalized Least-Mean Square Algorithm
%
% 	Estimator: xhat(n)=Y^{T}(n)thetahat(n-1)
%
%	thetahat is estimated using NLMS. 
%
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize xhat and thetahat

% Loop

for n=1:M,

	% Generate Y. Set elements of Y that does not exist to zero


	% Estimate of x



	% Update thetahat

	thetahat(n+1,:)=
end

% Shift thetahat one step so that row n corresponds to time n

thetahat=thetahat(2:M+1,:);
