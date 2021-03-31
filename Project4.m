%% Project 4 MA540
% Leah Rolf, Harley Hanes, James Savino
clear; clc; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
table=[58 59 60 61 62 63 64 65 66 67 68 69 70 71 72;115 117 120 123 126 129 132 135 139 142 146 150 154 159 164];
n=length(table);
covinit=[634.88 -235.04 21.66; -235.04 87.09 -8.03; 21.66 -8.03 0.74];
var=0.15;

weight=@(theta,height) theta(1)+theta(2)*(height/12)+theta(3)*(height/12).^2;
thetanom=[261.88;-88.18;11.96];
height=linspace(58,72,n);
epsilon=sqrt(var)*randn(1,n);
simulation=weight(thetanom,height)+epsilon;
doublesigma=2*var;

doublesigvec=[weight(thetanom,height)+doublesigma; weight(thetanom,height)-doublesigma];

% design matrix using Complex-Step
h=1e-16;

theta0com=[complex(thetanom(1),h) thetanom(2) thetanom(3)];
theta1com=[thetanom(1) complex(thetanom(2),h) thetanom(3)];
theta2com=[thetanom(1) thetanom(2) complex(thetanom(3),h)];

theta0=imag(weight(theta0com,height))/h;
theta1=imag(weight(theta1com,height))/h;
theta2=imag(weight(theta2com,height))/h;

design=[theta0;theta1;theta2]';

% Prediction interval for known interval
testxvalues=height(1:14)+.5;
interval=zeros(2,length(testxvalues));
for i=1:length(testxvalues)
heightextract=[1 (testxvalues(i)/12) (testxvalues(i)/12)^2]; 
predict=tcdf(0.975,n-3)*sqrt(var)*sqrt(1+heightextract*inv(design'*design)*heightextract');
yhatstar=heightextract*thetanom;
interval(:,i)=[yhatstar+predict;yhatstar-predict];
end

% Plots for known
figure('Renderer', 'painters', 'Position', [100 100 1050 500]);
subplot(1,2,1)
hold on
plot(height,weight(thetanom,height),'k','LineWidth',1)
scatter(table(1,:),simulation,'filled')
plot(height,doublesigvec,'r','LineWidth',1)
plot(testxvalues,interval,'b--','LineWidth',1)
axis([min(height) max(height) -inf inf])
xlabel('Height');ylabel('Weight')
legend({'Mean','Data','+2$\sigma$','-2$\sigma$','Predict'},'Interpreter','Latex','Location','Northwest')
title('Known values [58,72]')
hold off

subplot(1,2,2) % zooming in on same graph
hold on
plot(height,weight(thetanom,height),'k','LineWidth',1)
scatter(table(1,:),simulation,'filled')
plot(height,doublesigvec,'r','LineWidth',1)
plot(testxvalues,interval,'b--','LineWidth',1)
axis([height(4)-1 height(4)+1 simulation(4)-1 simulation(4)+1])
xlabel('Height');ylabel('Weight')
legend({'Mean','Data','+2$\sigma$','-2$\sigma$','Predict'},'Interpreter','Latex','Location','Northwest')
title('Known values [60,62]')
hold off
saveas(gcf,'Figures/Knowninterval.png')

% Unknown interval
heightpred=50:80;
undoublesig=[weight(thetanom,heightpred)+doublesigma;weight(thetanom,heightpred)-doublesigma];

% Prediction for unknown
predictxvalues=[50:79]+.5; % exploring 50-80 inches
intervalpred=zeros(2,length(predictxvalues));
for i=1:length(predictxvalues)
heightextract=[1 (predictxvalues(i)/12) (predictxvalues(i)/12)^2]; 
predict=tcdf(0.975,n-3)*sqrt(var)*sqrt(1+heightextract*inv(design'*design)*heightextract');
yhatstar=heightextract*thetanom;
intervalpred(:,i)=[yhatstar+predict;yhatstar-predict];
end

% Plots for unknown
figure('Renderer', 'painters', 'Position', [100 100 1050 500]);
subplot(1,2,1)
hold on
plot(heightpred,weight(thetanom,heightpred),'k','LineWidth',1)
scatter(height,simulation,'filled')
plot(heightpred,undoublesig,'r','LineWidth',1)
plot(predictxvalues,intervalpred,'b--','LineWidth',1)
axis([min(heightpred) max(heightpred) -inf inf])
xlabel('Height');ylabel('Weight')
legend({'Mean','Data','+2$\sigma$','-2$\sigma$','Predict'},'Interpreter','Latex','Location','Northwest')
title('Unknown values [50,80]')
hold off

subplot(1,2,2) % zooming in on same graph
hold on
plot(heightpred,weight(thetanom,heightpred),'k','LineWidth',1)
plot(heightpred,undoublesig,'r','LineWidth',1)
plot(predictxvalues,intervalpred,'b--','LineWidth',1)
axis([74 76 175 185])
xlabel('Height');ylabel('Weight')
legend({'Mean','+2$\sigma$','-2$\sigma$','Predict'},'Interpreter','Latex','Location','Northwest')
title('Unknown values [74,76]')
hold off
saveas(gcf,'Figures/Unknowninterval.png')

%% Problem 2

%% Problem 3

%% Problem 4
