function [xhat,num,den]=cw(y,PhixyNum,PhixyDen,PhiyyNum,PhiyyDen,m)

% CW	Causal Wiener filter
%
% [xhat,num,den]=cw(y,PhixyNum,PhixyDen,PhiyyNum,PhiyyDen,m)
%
% Input arguments:
%
% 	y 		observed variable
%
%	m		Prediction horizon, i.e. xhat(n+m)= h * y(n)
%
%	PhixyNum	Numerator of cross-spectrum Phixy
%	PhiuxyDen	Denominator of cross-spectrum Phixy
%	PhiyyNum	Numerator of spectrum Phiyy
%	PhiuyyDen	Denominator of spectrum Phiyy
%
% Output arguments
%
%	xhat		Estimate
%
%	num,den		Causal Wiener filter H=num/den 
%			The format is one-sided z-transform in ascending powers of 
%			z^{-1}, the first coefficient is a constant


% First deal with Phiyy

% Make numerator/denominators have equal length

[PhiyyNum,PhiyyDen]=eqsize(PhiyyNum,PhiyyDen);

% Do spectral factorization of Phiyy

PhiyyNumplus = specfac(PhiyyNum);
PhiyyDenplus = specfac(PhiyyDen);

% Create the negative part

PhiyyNumneg=polystar(PhiyyNumplus);
PhiyyDenneg=polystar(PhiyyDenplus);

% Add zeros so that the format is in ascending powers of z^{-1}

nb=length(PhiyyNumneg);
na=length(PhiyyDenneg);

PhiyyNumneg=[zeros(1,na-nb) PhiyyNumneg];
PhiyyDenneg=[zeros(1,nb-na) PhiyyDenneg];

% Now deal with Phixy

% Add z^{m} to Phixy

PhixyNum=[PhixyNum];
PhixyDen=[zeros(1,m) PhixyDen];

% Next make numerator/denominators have equal length

[PhixyNum,PhixyDen]=eqsize(PhixyNum,PhixyDen);

% Compute z^{m}Phi_{xy}/Phi_{yy}^{-1}

num=conv(PhixyNum,PhiyyDenneg);
den=conv(PhixyDen,PhiyyNumneg);

% Take positive part

[numpos,denpos]=pos(num,den);

% Multiply with inverse of Phi_{yy}^{+}

num=conv(PhiyyDenplus,numpos);
den=conv(PhiyyNumplus,denpos);

[num,den]=eqsize(num,den);
[num,den]=minreal(num,den);

% Filter

xhat=filter(num,den,y);

% Adjust so that xhat(k) is an estimate of x(k)

xhat=[zeros(m,1);xhat(1:length(xhat)-m)];

