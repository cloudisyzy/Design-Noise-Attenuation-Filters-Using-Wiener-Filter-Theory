function [numpos,denpos]=pos(num,den)

% POS	Extract the positive part of a rational function
%
% [numpos,denpos]=pos(num,den)
%
% num=[b_{0} b_{1} ... b_{nb}]
% den=[a_{0} a_{1} ... a_{na}]
%
% The function is given by 
%
%		(b_{0}+b_{1}z^{-1}+ ... b_{nb}z^{-nb})
%	H(z)=   --------------------------------------
%		(a_{0}+a_{1}z^{-1}+ ... a_{na}z^{-na})
% 
% Output arguments:
%
% numpos/denpos is the positive part of the system in the 
% same format as num/den


% convert from double-sided transform into one-sided transform (powers of z^{-1})
% + additional (anti)-delays

[num,den,dly]=delay(num,den);

if dly<0,

	% Convert to powers of z
	nnum=length(num); nden=length(den); nmax=max(nnum,nden);
	num=[num zeros(1,nmax-nnum-dly)];
	den=[den zeros(1,nmax-nden)];

	% Remove powers of z	
	[q,r]=deconv(num,den);
	num=r((length(r)-nmax+1):length(r))+q(length(q))*den;

else

	num=[zeros(1,dly) num];

end

% length(num)= length(den) so now we are back to powers of z^{-1} again
% Remove trailing zeros
num=rmtzeros(num);
den=rmtzeros(den);

[r,p,k]=residuez(num,den);
rnew=[];
pnew=[];
for l=1:length(r),

	if abs(p(l))<1,	rnew=[rnew;r(l)]; pnew=[pnew;p(l)]; end

end

[numpos,denpos]=residuez(rnew,pnew,k);
if isempty(denpos), denpos=1; end
if isempty(numpos), numpos=0; end


