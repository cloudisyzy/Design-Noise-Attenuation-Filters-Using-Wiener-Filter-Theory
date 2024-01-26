% SigmaYxhat = xcovhat(x, y, N)
%
%	y, x			- Data sequences
%	N			- Size of covariance matrix
%
%  xcovhat: Estimates SigmaYx=E[Y(n)x(n)]
%
%		where 
%
%	   	Y(n)=[y(n) y(n-1) ... y(n-N+1)]^{T}
%
%     
%     Author: Ziyue Yang
%     Date: 2024.01.24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SigmaYxhat = xcovhat(x, y, N)

%     The length of iteration is constrainted by the minimum length between
%     x, y, both column vectors
    M = min(length(x), length(y));

    SigmaYxhat = zeros(N, 1);

    for k = 1:N
        SigmaYxhat(k) = y(1:M-k+1)' * x(k:M) / M;
    end

end