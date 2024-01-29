function spec_comp(A,sigma2,Anoise,sigma2noise,numnc,dennc,numc,denc,thetahatfir);

%
%  spec_comp(A,sigma2,Anoise,sigma2noise,numnc,dennc,numc,denc,thetahatfir);
%	
%	A, sigma2 		- Estimated or true parameters of x
%	Anoise, sigma2noise	- Estimated or true parameters of v
% 	numnc, dennc		- Non-causal Wiener filter
% 	numc, denc		- Causal Wiener filter
% 	thetahatfir		- FIR Wiener filter
%	
% 		 
%	
% Plot bodediagrams for different estimation methods obtained from 
% est_add
%

% Ensure white background, so the black line is visible
whitebg(gcf,'white')

N=length(thetahatfir);

w=linspace(0,pi);
[mags,phases,ws]=dbode(1,A,1,w);
[magv,phasev,wv]=dbode(1,Anoise,1,w);
[magnc,phasenc,wnc]=dbode(numnc,dennc,1,w);
[magc,phasec,wc]=dbode(numc,denc,1,w);
[magfir,phasefir,wfir]=dbode(thetahatfir',1,1,w);
ws = ws/pi*0.5;
wv = wv/pi*0.5;
wnc = wnc/pi*0.5;
wc = wc/pi*0.5;
wfir = wfir/pi*0.5;
semilogy(ws,mags.^2*sigma2,'b',wv,magv.^2*sigma2noise,':k', ...
	 wnc,magnc.^2,'-.r',wc,magc.^2,'--g',wfir,magfir.^2,'-m', 'LineWidth',1.5)
legend('Speech','Noise','Non-Causal Filter','Causal Filter',['FIR(' int2str(N) ') Filter '])
title('Power Spectrum Density')
xlabel('Normalized Frequency \nu')
ylabel('Magnitude')

