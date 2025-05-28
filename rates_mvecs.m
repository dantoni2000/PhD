figure

plot([20:10:60],[20:10:60].^-0.5,'o-r')
hold on
plot([20:10:60],(2.*[10:10:50]).^-0.5 + exp(-5.*[10:10:50].^0.25),'o-y')
hold on
plot([20:10:60],(1.*[15:10:55]).^-0.5 + exp(-5.*[10:10:50].^0.25),'o-g')
hold on
plot([20:10:60],(5.*[10:10:50]).^-0.5 + exp(-2.*[10:10:50].^0.25),'o-m')
hold on
plot([20:10:60],(8.*[4:10:44]).^-0.5 + exp(-2.*[4:10:44].^0.25),'o-k')

title('rate of convergences for different values of $s$, $N$ and $m$', 'Interpreter','latex')
xlabel('matvecs')
ylabel('rates')
legend('Nystrom S', 'Nystrom 1Hutch 5Lanczos','Nystrom 2Hutch 5Lanczos', 'Nystrom 5Hutch 2Lanczos', 'Nystrom 8Hutch 2Lanczos')
