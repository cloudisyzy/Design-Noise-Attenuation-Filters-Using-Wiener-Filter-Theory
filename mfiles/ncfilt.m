function y=ncfilt(num,den,u);

% Non-causal filtering 
%
% y=ncfilt(num,den,u);
%
%
% num=[b_{0} b_{1} ... b_{nb}]
% den=[a_{0} a_{1} ... a_{na}]
%
% The filter is given by 
%
%		(b_{0}+b_{1}z^{-1}+ ... b_{nb}z^{-nb})
%	H(z)=   --------------------------------------
%		(a_{0}+a_{1}z^{-1}+ ... a_{na}z^{-na})


%
% The filter is factorized into stable causal
% and anti-causal factors before the filtering
%

K=length(u);
y=zeros(K,1);

% convert into num/den * delay  form

[num,den,dly]=delay(num,den);

% Factorize denominator into causal and anti-causal factors

[dencausal,denanticausal]=fact(den);

ndenanticausal=length(denanticausal)-1;

% The following factorization is used in the filtering
%
%  (1/dencausal)   (1/denanticausal)


% Do causal filtering
ycausal=filter(num,dencausal,u);

% Do anticausal filtering

y=acfilt(1,denanticausal,ycausal);

% Finally adjust dly

if dly>=0,

	y=[zeros(dly,1);y(1:(K-dly))];

else

	y=[y((-dly+1):K);zeros(-dly,1)];

end


