function [PhiNum,PhiDen]=filtspec(B,A,sigma2)

% FILTSPEC		Spectrum of filtered white noise
%
% [PhiNum,PhiDen]=filtspec(B,A,sigma2)
%
% Computes the spectrum of the signal
%
% 	y=B/A e
%
% where e is zero mean white noise with variance sigma2
%
% 
% B=[b_{0} b_{1} ... b_{nb}]
% A=[a_{0} a_{1} ... a_{na}]
%
% The filter is given by 
%
%		(b_{0}+b_{1}z^{-1}+ ... b_{nb}z^{-nb})
%	H(z)=   --------------------------------------
%		(a_{0}+a_{1}z^{-1}+ ... a_{na}z^{-na})
%
% PhiNum/PhiDen is the resulting spectrum in the same format
% as B/A
%


% Remove delays 

[B,A,dly]=delay(B,A);

nb=length(B)-1;
na=length(A)-1;

% Reverse order of coefficients

Ar=polystar(A);
Br=polystar(B);

% Convolution

PhiDen=conv(A,Ar);
PhiNum=conv(B,Br);


% Add zeros so that the format is in ascending powers of z^{-1}

PhiNum=[zeros(1,na-nb) PhiNum];
PhiDen=[zeros(1,nb-na) PhiDen];

PhiNum=sigma2*PhiNum;
