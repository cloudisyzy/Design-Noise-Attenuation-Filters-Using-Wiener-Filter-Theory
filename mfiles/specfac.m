function H = specfac(Phi)

% SPECFAC	Computes the spectral factorization
%
% H=specfac(Phi)
%
%  For a spectrum Phi(z)=sum{k=-n..n} a_{k} z^{-k}
%  the input argument should be Phi=[a_{-n} a_{-n+1} ... a_{n}]
%
%  The output H is the spectral factor H=[b_{0} b_{1} ... b_{n}]
%  such that H(z)=sum{k=0..n} b_{k}z^{-k} gives
%
%      Phi(z)=H(z)H(z^{-1})
%
%  The method is based on the calculation of the roots
%  and was developed by Peterka.
%
% The tolerance is 1e-10.
%

m=Phi;
if nargin == 1, 
    tol = 1e-10; 
    method = 0; 
end;

if nargin == 2,
    tol = 1e-10;
end;

rm = roots(m);
[n,bid] = size(rm);
rx  = [];

for i = 1:n
  if abs(rm(i))<1
    rx = [rx; rm(i)];
  end
end

[x,dum1] = zp2tf(rx,[],1);
xx = conv(x,polystar(x));


i = 0;
j = 1;
while i == 0,
  if xx(1,j) ~= 0,
    fact = m(1,j)/xx(1,j);
    x = sqrt(fact)*x;
    i = 1;
  end;
  j =j+1;
end;

if method == 1,
  while r>tol
        [x, rem]=deconv(m,x);
        lx=length(x);
        x=x(lx:-1:1)/x(lx);
        r=norm(rem,1)
  end;
  return;
end;

H=x;
