function [Num,Den]=add(Num1,Den1,Num2,Den2)

% ADD 	Add Z-transforms
%
% [Num,Den]=add(Num1,Den1,Num2,Den2)
%

% Make numerators and denominator transforms have equal length

[Num1,Den1]=eqsize(Num1,Den1);
[Num2,Den2]=eqsize(Num2,Den2);

% Add together

Num=conv(Num1,Den2)+conv(Den1,Num2);
Den=conv(Den1,Den2);

% Remove zeros

[Num,Den]=rmczeros(Num,Den);
