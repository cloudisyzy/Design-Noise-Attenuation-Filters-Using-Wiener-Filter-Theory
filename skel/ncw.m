% [xhat,num,den] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen)
%	
%	y			- y(n)=x(n)+v(n)
% 	PhixyNum,PhixyDen	- Cross-spectrum between x(n) and y(n)
% 	PhiyyNum,PhiyyDen	- Spectrum of y(n)
%	
% 	xhat		- Non-causal Wiener estimate of x(n) from y(n)
% 	num,den		- The Non-causal Wiener filter
%
%  ncw: Non-causal Wiener filtering.
%     
%     
%     Author: Ziyue Yang
%     Date: 2024.01.23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xhat,num,den] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen)

    addpath('./mfiles');

%     According to H(z)=Phi_xy (z)/ Phi_yy (z), see details in p.12 in
%     manual
    num = conv(PhixyNum, PhiyyDen);
    den = conv(PhixyDen, PhiyyNum);

%     ncfilt.m is used for non-causal filtering
    xhat = ncfilt(num, den, y);

end
