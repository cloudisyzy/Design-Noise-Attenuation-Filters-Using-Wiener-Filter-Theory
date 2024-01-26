function sigma2hat=rex1(lambdamuu,method)

% sigma2hat=rex1(lambdamuu,method)
%
%	lambdamuu	- Forgetting factor or step size
%	method		- 'rls', 'lms' or 'nlms' 
%	sigma2hat	- Normalized sum of the squared prediction error (y(n)-yhat(n|n-1)
%
%  rex1s	Generates plots for rex1
%
%		Requires the functions ar_id_lms, ar_id_nlms, ar_id_rls
%
%

clf
load rex1_all
if strcmp(method,'lms')

	[Ahat,yhat,sigma2hat]=ar_id_lms(y,2,lambdamuu);

elseif strcmp(method,'nlms')

	[Ahat,yhat,sigma2hat]=ar_id_nlms(y,2,lambdamuu);

elseif strcmp(method,'rls')

	[Ahat,yhat,sigma2hat]=ar_id_rls(y,2,lambdamuu);
%	sigma2hat=cov(y(5:length(y))-yhat(5:length(y)));

end

plot(Ahat(:,2))
hold on
plot(Ahat(:,3),'g--')
plot([1 K1],[A1(2) A1(2)])
plot([K1 K2-1],[A2(2) A2(2)])
plot([K2 K],[A3(2) A3(2)])
plot([1 K1],[A1(3) A1(3)],'g--')
plot([K1 K2-1],[A2(3) A2(3)],'g--')
plot([K2 K],[A3(3) A3(3)],'g--')
hold off
ylim=get(gca,'ylim');
set(gca,'ylim', [min(ylim(1), -1.5), ylim(2)])