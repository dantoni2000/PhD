% 08 04 2025 convergence of  Lanczos for the original and for the
% preconditioned matrix

clear all
close all
warning off

T = 1; dl = 20; s = 10; p = 0.01; mi = 0;

decadimento=7;

switch decadimento

    case 0
        n = 100;
        A = gen_gnp(n, 0.05);

    case 1
        n = 1000;
        A = zeros(n);
        for i=1:n-1
            A(i,i+1)=1; A(i+1,i)=1;
        end
        for i=1:50
            u = randi(n); v = randi(n);
            A(u,v)=1; A(v,u)=1;
        end


    case 2
        n = 1000;
        A = zeros(n);
        for i = 1:n-1
            A(i,i+1) = 1;
            A(i+1,i) = 1;
        end


    case 3
        n = 1000;
        A = zeros(n);
        for i = 1:n-1
            A(i,i+1) = 1;
            A(i+1,i) = 1;
        end
        A(1,end) = 1;
        A(end,1) = 1;

    case 4
        n = 1000;
        A = ones(n) - eye(n); % completo

    case 5
        n = 1000;
        A = zeros(n);
        hub = 1;
        A(hub, 2:1001) = 1; 
        A(2:1001, hub) = 1;

    case 6
        n = 1000;
        A = zeros(n);
        for i = 1:n-1
            A(i,i+1) = 1; A(i+1,i) = 1;
        end
        % aggiungo rami casuali
        p = 0.005;
        for i = 1:n
            for j = i+2:n
                if rand < p
                    A(i,j) = 1; A(j,i) = 1;
                end
            end
        end
        
    case 7
        n = 1000;
        k = 50; % cluster di grandi autovalori
        A = diag(ones(n,1)*1); % base identità
        A(1:k,1:k) = A(1:k,1:k) + 0.999*rand(k); % cluster vicino a 2
        A = (A + A')/2; % simmetrica

    case 8
        n = 1000;
        parts = [300,400,300];                % tre partizioni
        A = gen_random_multipartite(parts, 0.1); % probabilita p tra partizioni

end

[U,S,~] = svd(A);

% lambda_max = S(1,1);
% 
% S = diag(diag(S-S(1,1)));

% EE = sum(exp(lambda - lambda_max)) * exp(lambda_max);

f=@(x)(expm(x));

for t=1:T 
        
    for l=10:10:10*dl
        
        % Precondizionatore ottimale

        % L = [diag(S(1:l/2-1,1:l/2-1)); linspace(S(l/2-1,l/2-1),S(n-l/2+1,n-l/2+1),n-l+1)'; diag(S(n-l/2+1:end,n-l/2+1:end))];
        L = [zeros(n-l,1); diag(S(n-l+1:end,n-l+1:end)).^0.5];

        % L2 = [diag(S(1:l,1:l)).^0.5; zeros(n-l,1)];

        %Hutchinson sulla matrice precondizionata
        figure(1+l/10)
        [Pits(l/10),Ptr,h] = HUTCH_exp(A,mi,U,L,s,500,0,0);
    
        % hold on
        % [Pits2(l/10),Ptr2,h2] = HUTCH_exp(A,mi,U,L2,s,500,1,0);
    
        hold on
        [its(l/10),tr,k] = EST_function(A,mi,s,500,0,f);
        % f = semilogy(cA*rhoA.^[1:fix(its(l/10))+1],'m');

        legend([h h2 k],'$\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|$','roba','$\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|$','interpreter','latex','FontSize',12)
        %legend([h k],'$\frac{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|}{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v \right\|}$','$\frac{\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|}{\left\| v^T \log(A + \mu I) v \right\|}$','interpreter','latex','FontSize',12)
        xlabel('Lanczos Iterations')
        ylabel('Error')
    end

end