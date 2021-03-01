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


figure('Renderer', 'painters', 'Position', [100 100 800 600])
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

% Other design matrix using finite differences, for verification
% h=1e-6;
% psi1=(psi([alpha1+h,alpha11,alpha111],interval2)-psi(alpha0,interval2))/h;
% psi11=(psi([alpha1,alpha11+h,alpha111],interval2)-psi(alpha0,interval2))/h;
% psi111=(psi([alpha1,alpha11,alpha111+h],interval2)-psi(alpha0,interval2))/h;
% psidesign=[psi1; psi11; psi111]';
% design=psidesign;

% Solving for theta
thetaOLS=inv(design'*design)*design'*simulation';
psitheta=@(P) thetaOLS(1)*P.^2+thetaOLS(2)*P.^4+thetaOLS(3)*P.^6;
varianceOLS=(1/(n-3))*(simulation'-design*thetaOLS)'*(simulation'-design*thetaOLS);

% Formulating CI, not necessary but interesting
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
figure('Renderer', 'painters', 'Position', [100 100 1050 550])
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
legend({'Residuals','2$\sigma$'},'Interpreter','Latex')
hold off
saveas(gcf,'Figures/FittedandResiduals.png')

%% Problem 3

modelType='Generic SIR';
%Load SIR data
    D=load('SIR.txt');
    data.x=D(:,1);
    data.y=D(:,2);
    clear('D')
%Set SIR settings
paramsInit = [0.02; 0.15; 0.6];

%Get Optimal Paramestimates
[paramsOptimal,residualsOptimal,varEst,covEst] = LSQnonlin(...
    @(data,params)P3SIRoutput(data,params,'infected'),data,paramsInit);

%Plot Model Outputs over time
figure('Renderer', 'painters', 'Position', [100 100 1050 550])
subplot(1,2,1)
yOptimal=P3SIRoutput(data.x,paramsOptimal,'all');
plot(data.x,data.y,LineSpec(3,'marker'))
hold on
for i=1:3
    plot(data.x,yOptimal(:,i),LineSpec(i,'line'))
end
legend('data', 'Susceptible','Infected','Recovered','Interpreter','Latex')
xlabel('Time','Interpreter','Latex')
ylabel('Number of Individuals','Interpreter','Latex')
title(sprintf('SIR Model fit $(\\gamma=%.2g, \\delta=%.2g, r=%.2g)$',...
    paramsOptimal),'Interpreter','Latex')
subplot(1,2,2)
    hold on
    plot(data.x,residualsOptimal,LineSpec(1,'marker'))
    plot(data.x,ones(length(data.x),1)*2*sqrt(varEst),LineSpec(2,'line'))
    plot(data.x,-ones(length(data.x),1)*2*sqrt(varEst),LineSpec(2,'line'))
    legend('Residuals','$2\sigma$','Interpreter','Latex')
    axis([data.x(1) data.x(end) -max(abs(residualsOptimal)) max(abs(residualsOptimal))])



saveas(gcf,'Figures/SIRplot')


%Plot distributions and residuals
figure('Renderer', 'painters', 'Position', [100 100 1050 550])
%Get parameter distributions
paramSigma=sqrt(diag(covEst));
paramRanges=paramsOptimal+3*[-paramSigma paramSigma];
distPoints=100;
distRanges=NaN(3,100);
paramsPDF=NaN(3,100);
for i=1:3
    subplot(1,3,i)
    hold on
    distRanges(i,:)=linspace(paramRanges(i,1),paramRanges(i,2),distPoints);
    paramsPDF(i,:)=normpdf(distRanges(i,:),paramsOptimal(i),paramSigma(i));
    plot(distRanges(i,:),paramsPDF(i,:),LineSpec(i,'line'))
    if i==1
        ylabel('pdf')
        xlabel('$\gamma$','Interpreter','Latex')
    elseif i==2
        title(sprintf('pdf''s of parameter for %s model',modelType))
        xlabel('$\delta$','Interpreter','Latex')
    else
        xlabel('r','Interpreter','Latex')
    end  
end

