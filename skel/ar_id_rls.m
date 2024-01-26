function [Ahat,yhat,sigma2hat]=ar_id_rls(y,N,lambda)

% [Ahat,yhat,sigma2hat]=ar_id_rls(y,N,lambda)
%
%	y	  - Data sequence
%	N	  - Dimension of the parameter vector
%	lambda	  - Forgetting factor
%	Ahat	  - Matrix with estimates of the coefficients of the
%		    A polynomial [1, a1, ..., aN]. Row n corresponds 
%		    to the parameter estimates at time n.
%	yhat	  - One-step ahead prediction yhat(n)= yhat(n|n-1)
%	sigma2hat - Covariance of prediction error 
%
%
%  ar_id_rls: Recursive identification of an AR model using RLS
%
% 	
%     Author: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

