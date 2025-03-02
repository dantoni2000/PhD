%test 27 02 2025
%on mapping the s column of e.


err=1e-0;

C1 = zeros(20,1);
C1(1) = 1;
C1(2) = -1;
C2 = zeros(20,1);
C2(3) = 1;
C2(4) = -1;
C3 = zeros(20,1);
C3(1) = 1;
C3(2) = -1;
C3(3) = 1;
C3(4) = -1;
C = [C1 C2 C3];


E1 = zeros(20,1);
E1(1:2) = err*ones(2,1);

E2 = zeros(20,1);
E2(3:4) = err*ones(2,1);

E3 = zeros(20,1);
E3 = err*(E1+E2);

EJ = [E1 E2 E3];

ne = norm(EJ,'fro')^2;
nc = norm(C,'fro')^2;
comp = norm((eye(20,20)-(C+EJ)*pinv(C+EJ))*EJ,'fro')^2;
conf = ne*(nc^2 + ne*nc)/(nc^2 + 2*ne*nc +ne^2);

trueness = conf - comp;