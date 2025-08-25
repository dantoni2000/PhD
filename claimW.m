Omega = randn(n,l);
Omega1 = Omega(1:l,1:l);
Omega2 = Omega(l+1:end,1:l);

Z = Omega2*pinv(Omega1);
Inv = inv((eye(50)+sigma*0.5*(Z'*Z + (Z'*Z)')));

Zbrutta = diag(1./diag(almostG1(1:50,1:50))) *  Z;
InvBrutta = inv((eye(50)+sigma*0.5*(Zbrutta'*Zbrutta + (Zbrutta'*Zbrutta)')));

eigs(- InvBrutta + 10*Inv,50)