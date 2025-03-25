%test 15 03 2025 
%for a matrix C, we show that there exists Ej so that the bound is tight


    c1 = 1e-3; c2 = -0.5e-2; c3 = 4e-4; c4 = -1e-3;
    % c2 = c1; c3 = c1; c4 = c1;
    err = 0.005;

    C = [c1; c2; c3; c4];
    c = norm(C);
    theta = atan(-c1/c2);
    R = eye(4,4);
    R(1:2,1:2) = [cos(theta)  sin(theta);
                  -sin(theta) cos(theta)];
   
    RC = R*C;
    RE = [sqrt(err^2 -err^4/c^4*RC(2)^2 - err^4/c^4*RC(3)^2 -err^4/c^4*RC(4)^2); -err^2/c^2 * RC(2); -err^2/c^2 * RC(3); -err^2/c^2 * RC(4)];
    E = R'*RE;

    Ct = R'*(RC+RE);
    [Q,~] = qr(Ct,0);
    
    n1 = norm(Q*Q'*C)^2;
    n2 = norm(C)^2 - norm(E)^2;
    
    m = abs(n1 - n2);