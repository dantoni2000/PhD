% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 1; l = 50; S = 20;
MatVecsEST = zeros(S/5,1); MatVecsHutch = zeros(S/5,1); MatVecsinvHutch = zeros(S/5,1); 

decadimento=0;

switch decadimento
    case 0
        n = 2000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

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
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
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
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
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
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
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
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
    case 5
        n = 2000;
        mi = 1e-4;
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
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");
[its0,~] = EST(A,mi,S,500,1);

for t=1:T

    % Nystrom su A
    [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l);
    % P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
    P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
    trPREC = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));

    [iU,iLhat] = Nystrom(inv(A),l);  iLhat=inv(iLhat); ill = iLhat(l,l);
    iP = (ill+mi)^0.5*iU*(iLhat+mi*eye(l))^-0.5*iU' + (eye(n) - iU*iU');
    trPRECi = log((ill+mi)^(-l)*prod(diag(iLhat+mi*eye(l,l))));

       
    for s=5:5:S
        
        MatVecsHutch(s/5) = MatVecsHutch(s/5) + l;
        MatVecsinvHutch(s/5) = MatVecsinvHutch(s/5) + l;
        % trPREC(l/10,1) = trace(logm((P*P)^-1));
        % true_value(l/10,1) = trace(logm(P05*(A+mi*eye(n))*P05));
        
        %Hutchinson sulla matrice precondizionata
        [Pits,Ptr] = Preconditioned_HUTCH(A,mi,P,s,500,1,1);
        MatVecsHutch(s/5) = MatVecsHutch(s/5) + Pits*s;
        PreHutch(s/5,1) = trPREC + Ptr;
           
        %Hutchinson sulla matrice precondizionata (inversa)
        [iPits,iPtr] = Preconditioned_HUTCH(A,mi,iP,s,500,1,1);
        MatVecsinvHutch(s/5) = MatVecsinvHutch(s/5) + iPits*s;
        PreHutch(s/5,1) = trPRECi + iPtr;

        %stimatore direttamente su A
        z = fix(MatVecsHutch(s/5)/its0);
        [its,tr] = EST(A,mi,z,500,1);
        MatVecsEST(s/5) = its*z;
        HutchA(s/5,1) = tr;

         invz = fix(MatVecsinvHutch(s/5)/its0);
        [invits,invtr] = EST(A,mi,invz,500,1);
        MatVecsinvEST(s/5) = invits*invz;
        HutchinvA(s/5,1) = invtr;

        %Nystrom più grande: not worth in some cases
        BIGL = fix(MatVecsHutch(s/5));
        [BIGU,BIGLhat] = Nystrom(A,BIGL); LL = BIGLhat(BIGL,BIGL);
        trBIGL(s/5,1) = trace(logm(BIGU*BIGLhat*BIGU' + mi*eye(n)));

        % BIGL = 1900;
        % [BIGU,BIGLhat] = Nystrom(A,BIGL); LL = BIGLhat(BIGL,BIGL);
        % trBIGL(s/5,1) = trace(logm(BIGU*BIGLhat*BIGU' + mi*eye(n)));


        BIGinvL = fix(MatVecsinvHutch(s/5));
        [BIGinvU,BIGinvLhat] = Nystrom(inv(A),BIGinvL); BIGinvLhat = inv(BIGinvLhat); invLL = BIGinvLhat(BIGinvL,BIGinvL);
        trBIGinvL(s/5,1) = trace(logm(BIGinvU*BIGinvLhat*BIGinvU' + mi*eye(n)));

        % BIGinvL = 1900;
        % [BIGinvU,BIGinvLhat] = Nystrom(inv(A),BIGinvL); BIGinvLhat = inv(BIGinvLhat); invLL = BIGinvLhat(BIGinvL,BIGinvL);
        % trBIGinvL(s/5,1) = trace(logm(BIGinvU*BIGinvLhat*BIGinvU' + mi*eye(n)));

        
    end
    [sMatVecsEST,ord] = sort(MatVecsEST);


    figure(2)
    plot(sMatVecsEST,abs(PreHutch(ord)-trA)./abs(trA),'r')
    hold on
    plot(sMatVecsEST,abs(HutchA(ord)-trA)./abs(trA),'o-g')
    legend('PreHutch++','est(log(A+$\mu$ I))','trlog(A+$\mu$ I)','interpreter','latex')
    xlabel('matvecs')
    ylabel('error')
    title('error wrt matvecs by fixing Nystrom')
end