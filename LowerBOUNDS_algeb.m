% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off
T = 10;
n = 1000;
g = linspace(1,n,n)';
mi = 1;
ExperrtrBig=zeros(9,9);
Experrtr1=zeros(9,9);

for l=15:10:95

    for j = 1:9

        alpha = .5 +j;
        A = diag(g.^-alpha);
        trA(j) = sum(log(diag(A+mi*eye(n,n))),"all");
        mvecs((l-5)/10 ) = l;
        
        for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = Nystrom(A,l);
        ExperrtrBig((l-5)/10,j) = ExperrtrBig((l-5)/10,j) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA(j)).^2;
    
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        Experrtr1((l-5)/10,j) = Experrtr1((l-5)/10,j) + abs(trr - trA(j)).^2;
    
        for s = 4:2:l-6
            k = s;
            p = l - 5 - s;
    
            % boundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k+1) .* (sqrt(1-alpha.^(2*(n-k)))./sqrt(1-alpha.^2) + (k-1)/p .*sqrt(1-alpha.^(n-k))./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .*(1-alpha.^(n-k))./ (1-alpha));        
            LOWERboundPrecSTE(:,s/2-1) = (1+ (1+A(1,1)).^2 .*A(1,1) .* min(A(1,1),((1+2*k/(p-1)).*A(k+1,k+1)+(2*exp(1)^2*(k+p)/(p^2-1).* norm(A(k+1:n,k+1:n).^0.5,'fro')^2))) +2.*(1+A(1,1)).*min(A(1,1),((1+2*k/(p-1)).*A(k+1,k+1)+(2*exp(1)^2*(k+p)/(p^2-1).* norm(A(k+1:n,k+1:n).^0.5,'fro')^2)))).*(1+A(1,1)).^2;
            LOWERboundBIGNys(:,s/2-1) = ( (1+(1+A(1,1)).* min(A(1,1),((1+2*(k+5)/(p-1)).*A(k+6,k+6)+(2*exp(1)^2*(k+5+p)/(p^2-1).* norm(A(k+6:n,k+6:n).^0.5,'fro')^2 ))) ).*(1+A(1,1)) ).^2; 

        end

        %boundPrecSTE((l-5)/10,j) = 2 .* norm(A(l-4:n,l-4:n),'fro')^2;
        boundPrecSTE((l-5)/10,j) = (min(LOWERboundPrecSTE')).^-1 .* 2 .* norm(A(l-4:n,l-4:n),'fro')^2;
        boundBIGNys((l-5)/10,j) = (min(LOWERboundBIGNys')).^-1 .* norm(A(l+1:n,l+1:n).^0.5,'fro')^4;
        
        end
    
    end

end

for r=1:1:9
    
    for j = 1:8

    if (boundPrecSTE(r,j)' - boundBIGNys(r,j)')<=0 && (boundPrecSTE(r,j+1)' - boundBIGNys(r,j+1)')>=0
        nu(r) = 1.5 + (alpha-1.5) * ((j-1)/9);
        break    
    
    else 
        nu(r) = alpha;

    end
    
    end
    
end

etrBig = 1/T * ExperrtrBig;
etr1 = 1/T * Experrtr1;

for r=1:1:9
    figure(r)
    semilogy(1.5:1:9.5,etr1(r,:)','r')
    hold on
    semilogy(1.5:9.5,boundPrecSTE(r,:)','o-r')
    hold on
    semilogy(1.5:1:9.5,etrBig(r,:)','b')
    hold on
    semilogy(1.5:9.5,boundBIGNys(r,:)','o-b')

    xlabel('$\nu$', 'Interpreter','latex')
    ylabel('error')
    legend('error Nystrom+Lanczos', 'bound Nystrom+Lanczos', 'error Big Nystrom', 'bound Big Nystrom')

end

figure(100)
    % plot(10:10:100,nu,'-*r')
    % hold on
    % plot(10:10:100,badnu,'-*b')
    % hold on
    plot(nu,'-*k')
    xlabel('$k$', 'Interpreter', 'latex')
    ylabel('$\nu$', 'Interpreter', 'latex')
    title('Coefficient $\nu$ for which $trace \log((I + x_{ii}^{-\nu}))$ can be better approximated with PrecSTE', 'Interpreter','latex')
