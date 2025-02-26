%test 20 02 2025
%on best error Ej: the orthogonal one! 
%Hence, given C, properly choose two indeces of C and switch them and their sign!

% C = [ones(10,1); zeros(10,1)];
% 
% for j=10:2:20
% Ej = [zeros(20-j,1); ones(j,1)];
% 
% conf = norm(C)^2;
% comp((j-8)/2) = norm((eye(20,20)-(C+Ej)*pinv(C+Ej))*Ej)^2;
% end

err=1e-03;

C = [ones(10,1); -ones(10,1)];
Ej = err*ones(20,1);

ne = norm(Ej)^2;
nc = norm(C)^2;
comp = norm((eye(20,20)-(C+Ej)*pinv(C+Ej))*Ej)^2;
conf = ne*(nc^2 + ne*nc)/(nc^2 + 2*ne*nc +ne^2);

trueness = conf - comp;


C = [ones(10,1); -ones(10,1)];

E2j = err*zeros(20,1); E2j(10,1)=1; E2j(11,1)=1;
E2j = E2j*ne/norm(E2j);

nc = norm(C)^2;
comp2 = norm((eye(20,20)-(C+Ej)*pinv(C+Ej))*Ej)^2;
conf2 = ne*(nc^2 + ne*nc)/(nc^2 + 2*ne*nc +ne^2);
trueness2 = conf - comp;