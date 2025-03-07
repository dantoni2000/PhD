%test 27 02 2025
%on best error Ej: the orthogonal one! 

% C = [ones(10,1); zeros(10,1)];
% 
% for j=10:2:20
% Ej = [zeros(20-j,1); ones(j,1)];
% 
% conf = norm(C)^2;
% comp((j-8)/2) = norm((eye(20,20)-(C+Ej)*pinv(C+Ej))*Ej)^2;
% end

<<<<<<< HEAD
err=1e-0;
=======
err = 1e-04;
>>>>>>> 1d5be9a (initial commit)

C1 = zeros(20,1);
C1(1) = 1;
C1(2) = -1;

C2 = 2*C1;
C3 = 3*C1;
% C2 = zeros(20,1);
% C2(2) = 1;
% C2(3) = -1;
C = [C1 C2 C3];


E1 = zeros(20,1);
E1(5:12) = err*ones(8,1);

E2 = zeros(20,1);
E2(5:12) = 2*err*ones(8,1);

E3 = zeros(20,1);
E3(5:12) = 3*err*ones(8,1);
EJ = [E1 E2 E3];

ne = norm(EJ,'fro')^2;
nc = norm(C,'fro')^2;
comp = norm((eye(20,20)-(C+EJ)*pinv(C+EJ))*EJ,'fro')^2;
conf = ne*(nc^2 + ne*nc)/(nc^2 + 2*ne*nc +ne^2);

trueness = conf - comp;