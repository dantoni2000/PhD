clear all
close all
warning off

T = 1; dl = 20; s = 10;

decadimento=1;

switch decadimento
    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-2;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-2;
        G = diag(1./sqrt(g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(exp(-g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 4
        n = 1000;
        mi = 1e-4;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=diag(diag(G));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 5
        n = 1000;
        mi = 5e-6;
        alpha = 1; nu = 5/2;
        kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
        data_matrix = 1/n*randn(n,n);

        for row = 1:n

            for column = 1:n

                if row == column
                
                    A(row,column) = sqrt(pi)*gamma(nu)/(gamma(nu + 1/2)*alpha^(2*nu));

                else

                    A(row,column) = kernel(data_matrix(:,row),data_matrix(:,column));

                end

            end
        
        end

        [Q,G] = svd(A);
        G=diag(diag(G));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
end

% lgA = log(diag(G+mi*eye(n)));
lgA = logm(A+mi*eye(n));
ct = 0;
for l = [10 50 100 150 200]
    ct = ct+1;
    % Nystrom su inv(A)
    [iU,iLhat] = Nystrom(inv(A),l);  iLhat=inv(iLhat); ill = iLhat(l,l);
    iP = (ill+mi)^0.5*iU*(iLhat+mi*eye(l))^-0.5*iU' + (eye(n) - iU*iU');
    lgPiAPi = logm(iP*(A+mi*eye(n))*iP);

    % Nystrom su A
    [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l);
    P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
    lgPAP = logm(P*(A+mi*eye(n))*P);

    % combinazione dei due Nystrom
    V = [U(:,1:l/2) iU(:,1:l/2)];
    Lambdahat = diag([diag(Lhat(1:l/2,1:l/2)); diag(iLhat(1:l/2,1:l/2))]);
    Q = mi^0.5*V*(Lambdahat+mi*eye(l))^-0.5*V' + (eye(n) - V*V');
    lgQAQ = logm(Q*(A+mi*eye(n))*Q);

    N = 50000;

    for j=1:N
        v = randn(n,1);
        esA_G(j) = v'*lgA*v;
        esPiAPi_G(j) = v'*lgPiAPi*v;
        esPAP_G(j) = v'*lgPAP*v;
        esQAQ_G(j) = v'*lgQAQ*v;

        w = randsrc(n,1);
        esA_R(j) = w'*lgA*w;
        esPiAPi_R(j) = w'*lgPiAPi*w;
        esPAP_R(j) = w'*lgPAP*w;
        esQAQ_R(j) = w'*lgQAQ*w;

        z = unifrnd(-1,1,n,1);
        z = sqrt(n)*z/norm(z);
        esA_U(j) = z'*lgA*z;
        esPiAPi_U(j) = z'*lgPiAPi*z;
        esPAP_U(j) = z'*lgPAP*z;
        esQAQ_U(j) = z'*lgQAQ*z;
    end
    
    meanA_G = sum(esA_G(:))/N;
    varA_G(ct,1) = sum((esA_G-meanA_G).^2)/N;
    ctrlA_G(ct,1) = 2*norm(lgA,'fro')^2;
    
    meanPiAPi_G = sum(esPiAPi_G(:))/N;
    varPiAPi_G(ct,1) = sum((esPiAPi_G-meanPiAPi_G).^2)/N;
    ctrlPiAPi_G(ct,1) = 2*norm(lgPiAPi,'fro')^2;

    meanPAP_G = sum(esPAP_G(:))/N;
    varPAP_G(ct,1) = sum((esPAP_G-meanPAP_G).^2)/N;
    ctrlPAP_G(ct,1) = 2*norm(lgPAP,'fro')^2;

    meanQAQ_G = sum(esQAQ_G(:))/N;
    varQAQ_G(ct,1) = sum((esQAQ_G-meanQAQ_G).^2)/N;
    ctrlQAQ_G(ct,1) = 2*norm(lgQAQ,'fro')^2;
    

    meanA_R = sum(esA_R(:))/N;
    varA_R(ct,1) = sum((esA_R-meanA_R).^2)/N;
    ctrlA_R(ct,1) = 2*(norm(lgA,'fro')^2-norm(diag(lgA))^2);

    meanPiAPi_R = sum(esPiAPi_R(:))/N;
    varPiAPi_R(ct,1) = sum((esPiAPi_R-meanPiAPi_R).^2)/N;
    ctrlPiAPi_R(ct,1) = 2*(norm(lgPiAPi,'fro')^2-norm(diag(lgPiAPi))^2);

    meanPAP_R = sum(esPAP_R(:))/N;
    varPAP_R(ct,1) = sum((esPAP_R-meanPAP_R).^2)/N;
    ctrlPAP_R(ct,1) = 2*(norm(lgPAP,'fro')^2-norm(diag(lgPAP))^2);

    meanQAQ_R = sum(esQAQ_R(:))/N;
    varQAQ_R(ct,1) = sum((esQAQ_R-meanQAQ_R).^2)/N;
    ctrlQAQ_R(ct,1) = 2*(norm(lgQAQ,'fro')^2-norm(diag(lgQAQ))^2);


    meanA_U = sum(esA_U(:))/N;
    varA_U(ct,1) = sum((esA_U-meanA_U).^2)/N;
    ctrlA_U(ct,1) = (2*n/(2+n))*(norm(lgA,'fro')^2-1/n*trace(lgA)^2);

    meanPiAPi_U = sum(esPiAPi_U(:))/N;
    varPiAPi_U(ct,1) = sum((esPiAPi_U-meanPiAPi_U).^2)/N;
    ctrlPiAPi_U(ct,1) = (2*n/(2+n))*(norm(lgPiAPi,'fro')^2-1/n*trace(lgPiAPi)^2);

    meanPAP_U = sum(esPAP_U(:))/N;
    varPAP_U(ct,1) = sum((esPAP_U-meanPAP_U).^2)/N;
    ctrlPAP_U(ct,1) = (2*n/(2+n))*(norm(lgPAP,'fro')^2-1/n*trace(lgPAP)^2);

    meanQAQ_U = sum(esQAQ_U(:))/N;
    varQAQ_U(ct,1) = sum((esQAQ_U-meanQAQ_U).^2)/N;
    ctrlQAQ_U(ct,1) = (2*n/(2+n))*(norm(lgQAQ,'fro')^2-1/n*trace(lgQAQ)^2);

end

figure(2)
plot([10 50 100 150 200],varA_G,'o-g')
hold on
plot([10 50 100 150 200],ctrlA_G,'g')
hold on
plot([10 50 100 150 200],varPiAPi_G,'o-c')
hold on
plot([10 50 100 150 200],ctrlPiAPi_G,'c')
hold on
plot([10 50 100 150 200],varPAP_G,'o-r')
hold on
plot([10 50 100 150 200],ctrlPAP_G,'r')
hold on
plot([10 50 100 150 200],varQAQ_G,'o-b')
hold on
plot([10 50 100 150 200],ctrlQAQ_G,'b')
title('comparison of variances in case of Gaussian random vectors')
legend('emp var for $A$', 'true var for $A$', 'emp var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'true var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'emp var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'true var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'emp var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$', 'true var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$','interpreter','Latex','FontSize',10)
xlabel('rank of Nystrom')
ylabel('variance')

figure(3)
plot([10 50 100 150 200],varA_R,'o-g')
hold on
plot([10 50 100 150 200],ctrlA_R,'g')
hold on
plot([10 50 100 150 200],varPiAPi_R,'o-c')
hold on
plot([10 50 100 150 200],ctrlPiAPi_R,'c')
hold on
plot([10 50 100 150 200],varPAP_R,'o-r')
hold on
plot([10 50 100 150 200],ctrlPAP_R,'r')
hold on
plot([10 50 100 150 200],varQAQ_R,'o-b')
hold on
plot([10 50 100 150 200],ctrlQAQ_R,'b')

title('comparison of variances in case of Rademacher random vectors')
legend('emp var for $A$', 'true var for $A$', 'emp var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'true var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'emp var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'true var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'emp var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$', 'true var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$','interpreter','Latex','FontSize',10)
xlabel('rank of Nystrom')
ylabel('variance')

figure(4)
plot([10 50 100 150 200],varA_U,'o-g')
hold on
plot([10 50 100 150 200],ctrlA_U,'g')
hold on
plot([10 50 100 150 200],varPiAPi_U,'o-c')
hold on
plot([10 50 100 150 200],ctrlPiAPi_U,'c')
hold on
plot([10 50 100 150 200],varPAP_U,'o-r')
hold on
plot([10 50 100 150 200],ctrlPAP_U,'r')
hold on
plot([10 50 100 150 200],varQAQ_U,'o-b')
hold on
plot([10 50 100 150 200],ctrlQAQ_U,'b')
title('comparison of variances in case of Uniform random vectors')
legend('emp var for $A$', 'true var for $A$', 'emp var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'true var for $(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$', 'emp var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'true var for $(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$', 'emp var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$', 'true var for $(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$','interpreter','Latex','FontSize',10)
xlabel('rank of Nystrom')
ylabel('variance')