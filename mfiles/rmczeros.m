function [num2,den2]=rmczeros(num,den)

% RMCZEROS	Remove common zeros
%
% [num2,den2]=rmczeros(num,den)
%
% Example num=[0 0 4 5 6 0 0 0], den=[0 1 0 0 0 0]
% gives   num2=[0 4 5 6], den2=[1]

num=rmtzeros(num);
den=rmtzeros(den);

num2=num; den2=den;
if isempty(num2) | isempty(den2), return; end
if norm(num2)==0 & norm(den2)==0, num2=[]; den2=[]; return; end

while (num2(1)==0) & (den2(1)==0)
	num2=num2(2:length(num2)); 
	den2=den2(2:length(den2)); 
end

