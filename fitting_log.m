% x = linspace(0, 28, 10000);
% y = log(1+x);
% 
% % Fit polinomiale di grado 5
% p = polyfit(x, y, 9);
% y_fit = polyval(p, x);
% 
% figure;
% plot(x, y, 'b', 'LineWidth', 2); hold on;
% plot(x, y_fit, 'r--', 'LineWidth', 2);
% legend('Funzione originale', 'Polinomio fit');
% title('Approssimazione polinomiale (minimo quadrato, non minimax)');
% 
% err = max(abs(y_fit-y));



clc; clear; close all;

n = 9;  % grado polinomio
f = chebfun(@(x) log(1+x), [0 28]);

% --- Polinomio minimax di grado n ---
[p, E] = minimax(f, n);  % p è il polinomio, E è errore massimo

% --- Visualizzazione ---
xx = linspace(0,28,1000);
plot(xx, f(xx), 'b', 'LineWidth', 2); hold on;
plot(xx, p(xx), 'r--', 'LineWidth', 2);
legend('Funzione originale', 'Polinomio minimax');
title(['Polinomio minimax di grado ', num2str(n)]);

% --- Errore massimo ---
disp(['Errore massimo (minimax): ', num2str(E)]);








N = 1;
kA = 1; kA = 28;
cA = 2*sqrt(kA+1)*log(2*kA)*N;
rhoA = ((sqrt(kA+1)-1)/(sqrt(kA+1)+1))^2;

cA*rhoA^10