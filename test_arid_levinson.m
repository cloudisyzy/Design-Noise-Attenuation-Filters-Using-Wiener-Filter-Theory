clear

v = 0.4;
a = [1 0.1 -0.8 -0.27];
w = sqrt(v)*randn(15000,1);
x = filter(1,a,w);

[A, sigma2] = ar_id(x, 3)

[r,lg] = xcorr(x,'biased');
r(lg<0) = [];

[ar,e] = levinson(r,numel(a)-1)
