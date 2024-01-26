function y=acfilt(num,den,u)

% Anti-causal filtering (including a constant)
%
% y=acfilt(num,den,u)
%
%
% num=[b_{0} b_{1} ... b_{nb}]
% den=[a_{0} a_{1} ... a_{na}] 
% 
% are the numerator and denominator coeffiecients of the filter 
%
%		(b_{0}+b_{1}z^{-1}+ ... b_{nb}z^{-nb})    B(z)
%	H(z)=   -------------------------------------- = ------
%		(a_{0}+a_{1}z^{-1}+ ... a_{na}z^{-na})    A(z)
% 
% which is a stable anti-causal filter (possibly including a causal constant)
%
% Notice that the roots of the A polynomial must be outside the unit circle
% and that na>=nb for the filter to be anti-causal (including a constant) 

y=[];
num=rmtzeros(num);
den=rmtzeros(den);
if length(num)>length(den), error('Numerator order > Denominator order in anti-causal filtering') ; end

[num,den]=eqsize(num,den);

% Reverse filter coefficients
den=polystar(den);
num=polystar(num);

% Reverse input sequence

K=length(u);
u=u(K:-1:1);

% Filter

y=filter(num,den,u);

% Reverse filtered sequence

y=y(K:-1:1);

