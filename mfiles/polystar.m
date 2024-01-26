function pstar = polystar(p)

% POLYSTAR	Calculates p*(z)
%
% pstar = polystar(p)
%
%  p(z)=sum{k=0..n} a_{k}z^{-k}, i.e. p=[p_{0} ... p_{n}], gives
%
% p*(z)=sum{k=0..n} a_{n-k}z^{-k}, i.e. pstar=[p_{n} ... p_{0}]
%

[n,m] = size(p);
pstar = p(m:-1:1);
