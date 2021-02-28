%Project2 Main
clear; clc; close all
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1

%% Problem 2 
n=801;
alpha1=-389.4; alpha11=761.3; alpha111=61.5;
sigma=2.2; %sigma^2=4.84
 % Change for each case
interval2=linspace(0,.8,n);
epsilon=sigma*randn(1,n); 

psi=@(alpha,P) alpha(1)*P.^2+alpha(2)*P.^4+alpha(3)*P.^6;
P=.8; alpha0=[alpha1,alpha11,alpha111];
psifu=@(alpha) (psi(alpha0,P)+sigma*randn)-(alpha1*P.^2+alpha11*P.^4+alpha111*P.^6); 
[ALPHA,VALUE]=fminsearch(psifu,alpha0);

epsilon=sigma*randn(1,n); 
simulation=psi([alpha1 alpha11 alpha111],interval2)+randn(1,n)*sigma; % row vector

Rhat=simulation-psi(ALPHA,interval2);
varhat=(1/(n-3))*(Rhat)*(Rhat)';

figure(1)
hold on
scatter(interval2,simulation,15)
plot(interval2,psi(alpha0,interval2))
legend('Data w/ error','Fitted w/ nominal','Location','Northwest')
title('Nominal Model and Observed Values vs. Polarization')
xlabel('Polarization');ylabel('Helmholtz Energy')
hold off
saveas(gcf,'Figures/DataNominal.png')

% Part B
% Making design/data matrix
design1=interval2.^2'; design2=interval2.^4'; design3=interval2.^6';
design=[design1 design2 design3];

% Solving for theta
thetaOLS=inv(design'*design)*design'*simulation';
psitheta=@(P) thetaOLS(1)*P.^2+thetaOLS(2)*P.^4+thetaOLS(3)*P.^6;
varianceOLS=(1/(n-3))*(simulation'-design*thetaOLS)'*(simulation'-design*thetaOLS);

% Formulating CI, not necessary just interesting
xtxinv=inv(design'*design);
confidence2=zeros(3,2);

for i=1:3
    t1=thetaOLS(i);
    t2=tcdf(.975,158)*sigma*sqrt(xtxinv(i,i));
    confidence2(i,:)=[t1-t2,t1+t2];
end

psilower=@(P) confidence2(1,1)*P.^2+confidence2(2,1)*P.^4+confidence2(3,1)*P.^6+epsilon; 
psiupper=@(P) confidence2(1,2)*P.^2+confidence2(2,2)*P.^4+confidence2(3,2)*P.^6+epsilon;

% Plotting
figure(2)
subplot(1,2,1)
hold on
scatter(interval2,simulation)
plot(interval2,psitheta(interval2))
legend('Observed Data','Fitted Model','Location','Northwest')
xlabel('Polarization');ylabel('Helmholtz Energy')
axis([0 .8 -inf inf])
hold off

subplot(1,2,2)
hold on
scatter(interval2,simulation-psitheta(interval2),'filled')
plot(interval2, ones(1,n)*-2*sqrt(varianceOLS),'r')
plot(interval2, ones(1,n)*2*sqrt(varianceOLS),'r')
axis([0 .8 -inf inf])
legend({'Residuals','2$\hat{\sigma}$'},'Interpreter','Latex')
hold off
saveas(gcf,'Figures/FittedandResiduals.png')

%% Problem 3
%Load SIR data
    D=load('SIR.txt');
    data.x=D(:,1);
    data.y=D(:,2);
    clear('D')
