% [Ahat,sigma2hat]=ar_id(y,N)
%
%	y			- Data sequence
%	N			- Model order 
%	Ahat			- Estimate of AR polynomial [1, a1, ..., aN]
%	sigma2hat		- Estimate of noise variance 
%
%
%  ar_id: Identification of AR model
%
%         Model: y(n)+a_{1}y(n-1)+...+a_{N}y(n-N)=e(n)
%
%     
%     Author: Ziyue Yang
%     Date: 2024.01.25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ahat, sigma2hat] = ar_id(y, N)

%     theta = -R_YY^-1 * R_Yx, with x(n) = y(n+1).
    
    SigmaYY = covhat(y, N);
    SigmaYx = xcovhat(y(2:end), y, N);
    [xhat, theta] = firw(y, SigmaYx, SigmaYY);

%     Ahat = Ahat=[1, aˆ1, aˆ2, . . ., aˆN ], page98 in lecnotes
    Ahat = [1, -theta'];
    
%     y(n) , page.25, first N samples are abandoned since we use N history samples
%     to estimate a new sample
%     yhat(n|n-1) , page.25,
%     note that the index are misaligned. A easy way to illustrate:
%     y(3|2)    ->y(2|1)<-      y(1|0)      no y(0|-1), for 2nd order matlab index starts at 2
%     y(3)      ->y(2)<-        y(1)        y(0),       for 2nd order matlab index starts at 3
    y_n = y(N+1:end);
    yhat_n_n_minus_1 = xhat(N:end);

    M = min(length(y_n), length(yhat_n_n_minus_1));
    sigma2hat = sum((y_n(1:M) - yhat_n_n_minus_1(1:M)) .^ 2) / M;

end
