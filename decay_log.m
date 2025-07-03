close all

n = 100;

d = linspace(1.2,0.8,n);
g = linspace(1,n,n);
[Q,~] = qr(randn(n));

A = Q*diag(d)*Q';
Ag = Q*diag((0.8).^g)*Q';
Aa = Q*diag(g.^(-2))*Q';

subplot(2,2,1)
imagesc(log(abs(A+eye(n))), [-16,3])
colorbar

subplot(2,2,2)
imagesc(log(abs(logm(A+eye(n)))), [-16,3])
colorbar

l = 30;

[U,Lhat] = Nystrom(A,l);
Pl = U*(Lhat+eye(l))^-0.5*U' + (eye(n) - U*U');

subplot(2,2,3)
imagesc(log(abs(Pl*(A+eye(n))*Pl)), [-16,3])
colorbar

subplot(2,2,4)
imagesc(log(abs(logm(Pl*(A+eye(n))*Pl))), [-16,3])
colorbar


figure(2)

subplot(2,2,1)
imagesc(log(abs(Ag+eye(n))), [-16,3])
colorbar

subplot(2,2,2)
imagesc(log(abs(logm(Ag+eye(n)))), [-16,3])
colorbar

l = 30;

[Ug,Lghat] = Nystrom(Ag,l);
Pgl = Ug*(Lghat+eye(l))^-0.5*Ug' + (eye(n) - Ug*Ug');

subplot(2,2,3)
imagesc(log(abs(Pgl*(Ag+eye(n))*Pgl)), [-16,3])
colorbar

subplot(2,2,4)
imagesc(log(abs(logm(Pgl*(Ag+eye(n))*Pgl))), [-16,3])
colorbar


figure(3)

subplot(2,2,1)
imagesc(log(abs(Aa+eye(n))), [-16,3])
colorbar

subplot(2,2,2)
imagesc(log(abs(logm(Aa+eye(n)))), [-16,3])
colorbar

l = 30;

[Ua,Lahat] = Nystrom(Aa,l);
Pal = Ua*(Lahat+eye(l))^-0.5*Ua' + (eye(n) - Ua*Ua');

subplot(2,2,3)
imagesc(log(abs(Pal*(Aa+eye(n))*Pal)), [-16,3])
colorbar

subplot(2,2,4)
imagesc(log(abs(logm(Pal*(Aa+eye(n))*Pal))), [-16,3])
colorbar

figure
plot(eigs(A+eye(n),100))
hold on
plot(real(eigs(Pl*(A+eye(n))*Pl,100)))