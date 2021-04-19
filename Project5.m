%% Project 5
% Leah Rolf, Harley Hanes, James Savino
clear; clc; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
load q_data.txt;
q_data=sort(q_data);
M=15; k=8;
qm=2*rand(M,1);qm=sort(qm); 
qm=linspace(0,2,M);
latinhyp=2*lhsdesign(M,1); latinhyp=sort(latinhyp);
fun=@(q) (6*q.^2+3).*sin(6*q-4);

qregress=zeros(M,k); % solving for weights u
for i=1:k
    qregress(:,i)=q_data.^i;
end
X=[ones(M,1) qregress];
u=X\fun(q_data);

coeff=u;
poly=@(coeff,q) coeff(1)+coeff(2).*q+coeff(3).*q.^2+coeff(4.)*q.^3+coeff(5).*q.^4+coeff(6).*q.^5+coeff(7).*q.^6+coeff(8).*q.^7+coeff(9).*q.^8;
polynomial=poly(coeff,q_data); 

figure('Renderer', 'painters', 'Position', [100 100 750 500])
hold on
plot(qm, fun(qm))
plot(latinhyp, fun(latinhyp),'--')
plot(q_data,fun(q_data),':')
xlabel('q');ylabel('f(q)')
legend({'Uniform','Latin Hypercube','Q data'},'Location','Northwest')
title('Uniform, Latin Hypercube, Q data f(q)')
hold off
saveas(gcf,'Figures/Project5/1samples.png')

figure('Renderer', 'painters', 'Position', [100 100 750 500])
hold on
plot(q_data,fun(q_data))
plot(q_data,polynomial,'--')
xlabel('q');ylabel('f(q)')
legend({'Function','Surrogate'},'Location','Northwest')
title('f(q) and Surrogate f_s(q)') 
hold off
saveas(gcf,'Figures/Project5/1surrogate.png')

%%%%% Second interval
interval2=linspace(-.5,2.5,M);
polynomial2=poly(coeff,interval2);

figure('Renderer', 'painters', 'Position', [100 100 1000 500])
subplot(1,2,1)
hold on
plot(interval2,fun(interval2))
plot(interval2,polynomial2,'--')
axis([-.5 2.5 -inf inf])
xlabel('q');ylabel('f(q)')
legend({'Function','Surrogate'},'Location','Northeast')
title('f(q) and Surrogate f_s(q)') 
hold off

subplot(1,2,2)
hold on
plot(interval2,fun(interval2))
plot(interval2,polynomial2,'--')
axis([-.1 2.1 -16 26])
xlabel('q');ylabel('f(q)')
legend({'Function','Surrogate'},'Location','Southwest')
title('Functions from [-.1,2.1]') 
hold off

saveas(gcf,'Figures/Project5/1surrogate2.png')


%% Problem 2

%% Problem 3


