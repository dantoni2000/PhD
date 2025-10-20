% rho = 0.55;
% rho = 0.7;
% rho = 0.85;

clear all 

n = 1000; c = 1e+3; delta = 10^1; m = 5;

beta = linspace(0.95,0.55,9);
hatn = log(delta/c)./log(beta);

hatn = floor(hatn);


for j = 1:9

    for l = 5:1:n

        LHS(j) = 2*(m+1)*(log(1+c*beta(j).^(l-m+1)))^2;
        
        if hatn(j)>l
        
            RHS1(j) = ((hatn(j)-l-1)^2)*(log(1+delta))^2;
            RHS2(j) = c^2*(beta(j).^(2*(hatn(j)+1))./(1-beta(j)^2) + 2*beta(j).^(2*(hatn(j)+2))./((1-beta(j)^2)^2));
        
            RHS(j) = RHS1(j) + 1/(1+delta)^2*RHS2(j);
        
        else 
        
            RHS(j) = c^2*(beta(j).^(2*(l+2))./(1-beta(j)^2) + 2*beta(j).^(2*(l+3))./((1-beta(j)^2)^2));
        
        end
    
        if LHS(j)<=RHS(j)
            wl(j) = l;
            break 
        end

    end

end


wl(wl>120)=1000;
plot(beta,wl,'-b', 'LineWidth', 5)
xlabel('$\beta$','fontsize',18, 'Interpreter','latex')
ylabel('$l$','fontsize',18, 'Interpreter','latex')
title('Number of $l$ to beat strategy B','fontsize',18,'Interpreter','latex')