function [yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0)

%
% [yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0)
%	
%	y	- measured signal, y=[y(1)'; y(2)'; ...; y(K)']
%	F,G,H	- State space model
%	R1	- covariance of the process noise w(n)
%	R2	- covariance of the measurement noise v(n)
%	x0	- Initial state value
%	Q0	- Initial error covariance
%	
% 	yhat	- estimate of y(n+1)
% 	xhatfilt - filtered x estimate x(n|n)
% 	xhatpred - predicted x estimate x(n+1|n)
% 	P	- Final posterior error covariance
% 	Q	- Final error covariance
%	
%
%  kalman: Kalman filter for the system
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

[N,ny]=size(y);   % ny=dimension of y(n), N=length of the data sequence
[nx,nx]=size(F);  % nx=dimension of x(n)

xhatfilt=zeros(N,nx); % Preallocation
xhatpred=zeros(N,nx); % Preallocation
yhat=zeros(N,1);      % Preallocation

% Initial round. The initial round is made separately since
% all MATLAB arrays start at index one.

xhatpred(1,:)=(F*x0)';
P=F*Q0*F'+G*R1*G';
L=P*H'/(H*P*H'+R2);
xhatfilt(1,:)=xhatpred(1,:)+(y(1,:)-xhatpred(1,:)*H')*L';
Q=P-P*H'/(H*P*H'+R2)*H*P;
yhat(1,:)=xhatpred(1,:)*H'; % Better to calculate yhat in a single 
			    % vector operation at the end?

for n=1:N-1
  % Calculate xhat(n+1|n), xhat(n+1|n+1) and yhat(n+1|n) here!
  
  
end

