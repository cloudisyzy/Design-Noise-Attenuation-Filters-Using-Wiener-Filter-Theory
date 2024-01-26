function [x,y] = gen_ss(F, G, H, R1, R2, x0, K)

%
% [x,y] = gen_ss(F, G, H, R1, R2, x0)
%	
%	F,G,H	- State space model
%	R1	- covariance of the process noise w(n)
%	R2	- covariance of the measurement noise v(n)
%	x0	- Initial state value
%	K	- Number of points
%	
% 	y	- generated signal, y=[y(1)'; y(2)'; ...; y(K)']
% 	x	- generated states, x=[x(1)'; x(2)'; ...; x(K)']
%	
%
%  gen_ss:
%     Generate signal and state sequences for the system
%             x(n+1) = F x(n) + G w(n)
%               y(n) = H x(n) + v(n)
%     
%  All data is stored as a column of N rows where each row
%  is the transposed state/data vector at the corresponding
%  time instance.
%               
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
