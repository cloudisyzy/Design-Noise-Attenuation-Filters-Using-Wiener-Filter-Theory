% [xhat,theta] = firw(y, SigmaYx, SigmaYY)
%	
%	y	- y(n)=x(n)+v(n)
% 	SigmaYY	- E[Y(n) (Y(n))']
%	SigmaYx	- E[Y(n) x(n)]
%	
% 	xhat	- FIR Wiener estimate of x(n) from y(n)
% 	theta	- FIR Wiener filter.
%	
%
%  firw: FIR Wiener estimate of x(n) from y(n)
%     
%     
%     Author: Ziyue Yang
%     Date: 2024.01.23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xhat,theta] = firw(y, SigmaYx, SigmaYY)

    addpath('./mfiles');

%     The 2nd line executes faster and more accurate that 1st line.
    % theta = inv(SigmaYY) * SigmaYx;
    theta = SigmaYY\SigmaYx;  % SigmaYY--N*N; SigmaYx--N*1

%%    First realization of computing xhat using matrix operations, hard to implement but helps to understand FIR Wiener
%     get toeplitz form of y, e.g., y = [1,2,3,4]; y_toeplitz = 
%     [ 1     2     3     4
%       2     1     2     3
%       3     2     1     2
%       4     3     2     1 ]
%     Gamma, in page.10, which is aka Y(n)^T, the first N lines of Gamma, 
%     where N=max_k for covariance R(k)
%     should be the lower 
%     triangular part of y_toeplitz,
%     with upper part equals to 0 for all elements.
%     e.g., Gamma_first_N_lines = [ 1     0     0     0
%                                   2     1     0     0
%                                   3     2     1     0
%                                   4     3     2     1 ]
%     If the length of y is L, for the rest lines of Gamma should be 
%     the flipped vector of y(i-N+1:i), e.g., [y(10),y(9),y(8),y(7)] for
%     is the 10th (i=10) line of Gamma if N=4.

    L = length(y);
    N = length(SigmaYx); % Also, [N,N] = size(SigmaYY)
    
    y_toeplitz = toeplitz(y(1:N));
    Gamma_first_N_lines = tril(y_toeplitz);
    Gamma = zeros(L, N);

    for i = 1:L
        if i <= N
            Gamma(i,:) = Gamma_first_N_lines(i,:);
        else
            Gamma(i,:) = flip(y(i-N+1:i));
        end
    end

%     By defination, xˆ(n) = Gamma * θ
    xhat = Gamma * theta;
    
%%     Alternatively, use 'conv' or 'filter', which is way much easier:
%     xhat = conv(theta, y); %% or xhat = filter(theta, y, 'full')
%     xhat = xhat(1:length(y));
    
end