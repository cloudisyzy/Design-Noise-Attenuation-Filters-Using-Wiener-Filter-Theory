function [num2,den2,dly]=delay(num,den);

% DLY	Delay in transfer function
%
% [num,den,dly]=delay(num,den);
%
% Counts and removes time-shifts in transfer functions
%
% Example:  	num=[0 1], den=[1 2] gives num2=1; den2=[1 2]; dly=1;
%		num=[0 1], den=[0 0 1 2] gives num=1; den2=[1 2]; dly=-1;
%

num=num(:).'; den=den(:).'; % Ensure row vectors
dly=[]; num2=num; den2=den;
[num2,den2]=rmczeros(num,den);
if isempty(num2) | isempty(den2), return; end
dly=0;
if num2(1),

	while den2(1)==0, den2=den2(2:length(den2)); dly=dly-1; end

else

	while num2(1)==0, num2=num2(2:length(num2)); dly=dly+1; end

end

