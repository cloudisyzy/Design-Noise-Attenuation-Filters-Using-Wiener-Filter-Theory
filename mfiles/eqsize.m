function [x,y]=eqsize(u,w)

% EQSIZE	Equal size (one-sided transform version)
%
% [x,y]=eqlength(u,w)
%
% Make two vectors have equal size by extendening with
% trailing zeros. 

% Remove trailing zeros from u and w

u=rmtzeros(u);
w=rmtzeros(w);

% Determine length of filter 

nu=length(u);
nw=length(w);
n=max(nw,nu);

% Make the vectors have equal length

x=[u zeros(1,n-nu)];
y=[w zeros(1,n-nw)];

