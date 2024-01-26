function y=rmtzeros(u)

% RMTZEROS	Remove trailing zeros
%
% y=rmtzeros(u)
%
% Example: 	u=[0 1 3 0 0] gives y=[0 1 3]
%

y=u;
if isempty(y), return; end
if norm(y)==0, y=[]; return;end

while y(length(y))==0, y=y(1:length(y)-1); end
