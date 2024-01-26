function [Bstab,Bunstab]=fact(A)

% [Bstab,Bunstab]=fact(A)
%
% factorizes A into its stable and unstable parts
% roots on the unit circle are considered as unstable.
%
% The unstable part is monic, i.e it begins with a 1.

if length(A)==1,Bstab=A; Bunstab=1; return,end
A1=A(1); A=A/A1;
v = roots(A);

vus = find(abs(v)>=(1-10*eps));
vunstab=v(vus);
vs = find(abs(v)<(1-10*eps));
vstab=v(vs);
Bstab=real(A1*poly(vstab));
Bunstab=real(poly(vunstab));
