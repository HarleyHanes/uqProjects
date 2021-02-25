%Project2 Main
clear; clc;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1

%% Problem 2
alpha1=-389.4; alpha11=761.3; alpha111=61.5;
sigma=2.2; %sigma^2=4.84
n=161; % Change for each case
interval2=linspace(0,.8,n);

psi=@(P) alpha1*P.^2+alpha11*P.^4+alpha111*P.^6;
%[x1,x2,x3]=lsqnonlin(psi,[alpha1,alpha11,alpha111]);
% X=fmincon(@(P) x(1)*P.^2+x(2)*P.^4+x(3)*P.^6,[alpha1;alpha11;alpha111])

epsilon=sigma*randn(1,n); 
simulation=psi(interval2)+epsilon;

% Part B
% Making design/data matrix
design1=interval2.^2'; design2=interval2.^4'; design3=interval2.^6';
design=[design1 design2 design3];

% Solving for theta
thetaOLS=inv(design'*design)*design'*simulation';
psitheta=@(P) thetaOLS(1)*P.^2+thetaOLS(2)*P.^4+thetaOLS(3)*P.^6;
varianceOLS=(1/(n-3))*(simulation'-design*thetaOLS)'*(simulation'-design*thetaOLS);

% Formulating CI
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
hold off

subplot(1,2,2)
hold on
scatter(interval2,simulation-psitheta(interval2),'filled')
plot(interval2, ones(1,n)*-2*sigma,'r')
plot(interval2, ones(1,n)*2*sigma,'r')
axis([0 .8 -inf inf])
legend({'Residuals','2$\sigma$'},'Interpreter','Latex')
hold off

%% Problem 3
%Load SIR data
    D=load('SIR.txt');
    data.x=D(:,1);
    data.y=D(:,2);
    clear('D')
