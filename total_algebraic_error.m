% rho = 0.55;
% rho = 0.7;
% rho = 0.85;

clear all 

n = 1000; c = 1e+4; delta = 1; m = 5;

alpha = -linspace(2,8,13);

hatn = (delta/c).^(1./alpha);

hatn = floor(hatn);


for j = 1:13

    for l = 5:1:n

        LHS(j) = 2*(m+1)*(log(1+c*(l-m+1).^alpha(j)))^2;
        
        if hatn(j)>l
        
            RHS1(j) = ((hatn(j)-l-1)^2)*(log(1+delta))^2;
            RHS2(j) = c^2*((hatn(j)+2).^(2*alpha(j)+2)./((2*alpha(j)+1).*(alpha(j)+1)) + 1/(2*alpha(j)+1)*((2-1/(alpha(j)+1))*(n+1)^(2*alpha(j)+2)));
        
            RHS(j) = RHS1(j) + 1/(1+delta)^2*RHS2(j);
        
        else 
        
            RHS(j) = c^2*((l+3).^(2*alpha(j)+2)./((2*alpha(j)+1).*(alpha(j)+1)) + 1/(2*alpha(j)+1)*((2-1/(alpha(j)+1))*(n+1)^(2*alpha(j)+2)) );
        
        end
    
        if LHS(j)<=1.5*RHS(j)
            wl(j) = l;
            break 
        end

    end

end
    
    
plot(alpha,wl,'-b', 'LineWidth', 5)
xlabel('$\alpha$','fontsize',18, 'Interpreter','latex')
ylabel('$l$','fontsize',18, 'Interpreter','latex')
title('Number of $l$ to beat strategy B','fontsize',18,'Interpreter','latex')