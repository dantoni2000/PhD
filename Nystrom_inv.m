function [U,Lhat] = Nystrom_inv(A,l) %doesn't work...and leverage is unfeasible for large scale

[~,n] = size(A);
[~,L] = Nystrom(A,4);

shift = 1.5*L(1,1);
shA = shift*eye(n) - A;

[U,shiftLhat] = Nystrom(shA,l-4);

Lhat = -shiftLhat + shift*eye(l-4);

