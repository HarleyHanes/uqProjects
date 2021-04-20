%% Project 4 MA540
% Leah Rolf, Harley Hanes, James Savino
clear; clc; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
table=[58 59 60 61 62 63 64 65 66 67 68 69 70 71 72;115 117 120 123 126 129 132 135 139 142 146 150 154 159 164];
n=length(table);
V=[634.88 -235.04 21.66; -235.04 87.09 -8.03; 21.66 -8.03 0.74];
var=0.15;
V_sig = sqrt(var)*eye(n,n);

weight=@(theta,height) theta(1)+theta(2)*(height/12)+theta(3)*(height/12).^2;
thetanom=[261.88;-88.18;11.96];
height=linspace(58,72,n);
epsilon=var*randn(1,n);
simulation=weight(thetanom,height)+epsilon;

% design matrix using Complex-Step and analytical
h=1e-16;

theta0com=[complex(thetanom(1),h) thetanom(2) thetanom(3)];
theta1com=[thetanom(1) complex(thetanom(2),h) thetanom(3)];
theta2com=[thetanom(1) thetanom(2) complex(thetanom(3),h)];

theta0=imag(weight(theta0com,height))/h;
theta1=imag(weight(theta1com,height))/h;
theta2=imag(weight(theta2com,height))/h;

% X=[theta0;theta1;theta2]';
X=[ones(n,1) (height/12)' ((height/12).^2)'];

% mean, covariance, and sd
mean = X*thetanom;
R = weight(thetanom,height)' - X*thetanom;
S2 = (1/(n-3))*(R'*R);
S = sqrt(S2);
V = S2*inv(X'*X);
covY = X*V*X' + V_sig;
sd_bd = 2*sqrt(diag(covY));

doublesigvec=[mean+sd_bd mean-sd_bd];

% Prediction interval for known interval
testxvalues=height(1:14)+.5;
testX=[ones(length(testxvalues),1) (testxvalues/12)' ((testxvalues/12).^2)'];
interval=zeros(2,length(testxvalues));

for i=1:length(testxvalues)
    heightextract=[1 (testxvalues(i)/12) (testxvalues(i)/12)^2]; 
    predict=tinv(1-.05/2,length(testxvalues)-3)*sqrt(var)*sqrt(1+testX(i,:)*inv(X'*X)*testX(i,:)');
    yhatstar=weight(thetanom,testxvalues(i));
    interval(:,i)=[yhatstar+predict;yhatstar-predict];
end

% Plots for known
figure('Renderer', 'painters', 'Position', [100 100 1050 500]);
hold on
plot(height,mean,'k','LineWidth',1)
scatter(table(1,:),simulation,'filled')
plot(height,doublesigvec,'r','LineWidth',2)
plot(testxvalues,interval,'b--','LineWidth',2)
axis([min(height) max(height) -inf inf])
xlabel('Height (inches)');ylabel('Weight (lbs)')
legend({'Mean','Data','+2$\sigma_Y$','-2$\sigma_Y$','Prediction Interval'},'Interpreter','Latex','Location','Northwest')
title('Known values [58,72]')
hold off
% saveas(gcf,'Figures/Knowninterval.png')

% Unknown interval
heightpred=50:80;
undoublesig=[weight(thetanom,heightpred)+sd_bd(1); weight(thetanom,heightpred)-sd_bd(1)];

% Prediction for unknown
predictxvalues=[50:79]+.5; % exploring 50-80 inches
intervalpred=zeros(2,length(predictxvalues));
predX=[ones(length(predictxvalues),1) (predictxvalues/12)' ((predictxvalues/12).^2)'];
for i=1:length(predictxvalues)
    heightextract=[1 (predictxvalues(i)/12) (predictxvalues(i)/12)^2]; 
    predict=tinv(1-.05/2,length(predictxvalues)-3)*sqrt(var)*sqrt(1+predX(i,:)*inv(X'*X)*predX(i,:)');
    yhatstar=weight(thetanom,predictxvalues(i));
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
xlabel('Height (inches)');ylabel('Weight (lbs)')
legend({'Mean','Data','+2$\sigma_Y$','-2$\sigma_Y$','Prediction Interval'},'Interpreter','Latex','Location','Northwest')
title('Unknown values [50,80]')
hold off

subplot(1,2,2) % zooming in on same graph
hold on
plot(heightpred,weight(thetanom,heightpred),'k','LineWidth',1)
plot(heightpred,undoublesig,'r','LineWidth',1)
plot(predictxvalues,intervalpred,'b--','LineWidth',1)
axis([75 79 180 198])
xlabel('Height (inches)');ylabel('Weight (lbs)')
% legend({'Mean','+2$\sigma$','Predict','-2$\sigma$'},'Interpreter','Latex','Location','Northwest')
title('Unknown values [75,79]')
hold off
% saveas(gcf,'Figures/Unknowninterval.png')

%% Problem 2

%% Problem 3/4
clear
    %Initialize problem
        %Load SIR data
        D=load('SIR.txt');
        data.xdata=D(:,1);
        data.ydata=D(:,2);
        clear('D')
        %Set initial conditions for SIR
        inits=[1000-127.1233; 127.1233; 0];
        
        %State parameter names and ranges
        params={
                {'gamma',0.0091172 ,0,1}
                {'delta', 0.1941,0,1}
                {'r', 0.78621,0,1}
            };
        %Define error function
        model.ssfun=@(params,data)sirSS(data,params,inits,'no k');
        %Define parameter variances
        options.qcov=[3.996604557707984e-08,7.294589132883216e-08,2.080488101816833e-07;7.294589132883216e-08,9.379222001817734e-05,9.841326827550314e-05;2.080488101816833e-07,9.841326827550314e-05,1.720799701083470e-04];

        %Set observation variance
        model.sigma2=3.890110168320449e+02;
        %Set number of simulations
        options.nsimu=10000;
        options.updatesigma=1;
        nSample=1000;
        
    %Run DRAM
        [results,chain,s2chain]=mcmcrun(model,data,params,options);
    %Print Chain results
        cov(chain)
        chainstats(chain,results)
        fprintf('\nsigma^2 mean=%.4e\n',mean(s2chain))
        fprintf('V_{DRAM}-V_{OLS}')
        results.qcov-options.qcov
        fprintf('(V_{DRAM}-V_{OLS})./V_{DRAM}')
        (results.qcov-options.qcov)./results.qcov
    %Get distribution
      [bandwidth_g,density_g,gMesh,cdf_q]=kde(chain(:,1));
      [bandwidth_d,density_d,dMesh,cdf_d]=kde(chain(:,2));
      [bandwidth_r,density_r,rMesh,cdf_r]=kde(chain(:,3));
    %Get Prediction Intervals
        predResults=mcmcpred(results,chain,s2chain,data.xdata,...
            @(data,params)SIReval(data,params,inits),nSample);
        dramInt=[predResults.obslims{1}{1}(3,:);
            predResults.predlims{1}{1}(3,:);
            predResults.predlims{1}{1}(2,:);
            predResults.predlims{1}{1}(1,:);
            predResults.obslims{1}{1}(1,:)];
        plotIntervals(data,dramInt);
legend('$95\%$ Prediction','$95\%$ Credible','Model Result','Data','Interpreter','Latex')
        %title('DRAM Intervals')
        xlabel('Time')
        ylabel('Infected Individuals')
saveas(gcf,'Figures/P4_DRAMInterval.png')
%% Problem 4
params=[.0091172  0.1941 0.78621];
Sens=getJacobian(@(params)SIReval(data.xdata,params,inits),...
    params);
sigma2=426.8;
V=sigma2*inv(Sens'*(Sens));
predVar=Sens*V*Sens';
confVar=predVar+sigma2*eye(size(predVar));

baseInfec=SIReval(data.xdata,params,inits);
predInfecInt=baseInfec+([2 -2].*sqrt(diag(predVar)));
confInfecInt=baseInfec+([2 -2].*sqrt(diag(confVar)));
plotMat=[confInfecInt(:,1) predInfecInt(:,1) baseInfec predInfecInt(:,2) confInfecInt(:,2)]';
plotIntervals(data,plotMat)
legend('$2\sigma$ Prediction','$2\sigma$ Credible','Model Result','Data','Interpreter','Latex')
        %title('Linearized Model Intervals')
        xlabel('Time')
        ylabel('Infected Individuals')
saveas(gcf,'Figures/P4_LinearizedInterval.png')
        
%% Comparing 3/4
figure('Renderer', 'painters', 'Position', [100 100 600 500])
        
plot(data.xdata,dramInt(1,:)-plotMat(1,:),LineSpec(1,'marker'))
hold on
plot(data.xdata,dramInt(2,:)-plotMat(2,:),LineSpec(2,'marker'))
	xlabel('Time')
    ylabel('Linearized Interval Error')
    legend('Prediction Interval','Credible Interval')
saveas(gcf,'Figures/P4_LinearizedIntervalError.png')

%% Support Functions 

function ss=sirSS(data,params,inits,useK)
if strcmpi(useK,'no k')
    params=[params(1) 1 params(2) params(3)];
end
ode_options = odeset('RelTol',1e-6);
[~,Y]=ode45(@SIR_rhs,data.xdata,inits,ode_options,params);
ss=sum((data.ydata-Y(:,2)).^2);
end

function I=SIReval(t,params,inits)
ode_options = odeset('RelTol',1e-6);
params=[params(1) 1 params(2) params(3)];
[~,Y]=ode45(@SIR_rhs,t,inits,ode_options,params);
I=Y(:,2);
end

function plotIntervals(data,intervals)
time=data.xdata;
figure('Renderer', 'painters', 'Position', [100 100 600 500])
hold on
%Plot Confidence Intervals
fillyy(time',intervals(1,:)',intervals (5,:)','g')
fillyy(time',intervals(2,:)',intervals (4,:)','r')
% x = [time(1) time' time(end) fliplr(time)'];
% inBetween = [intervals(5,1) intervals(1,:) intervals(5,end)  intervals(end,:)];
% fill(x, inBetween,'g','Linestyle','none')
% Plot Prediction Interval
% inBetween = [intervals(4,1) intervals(2,:) intervals(4,end) intervals(4,:)];
% fill(x, inBetween,'r','Linestyle','none')
%Plot Estimate
plot(time,intervals(3,:),'LineWidth',1.8)
plot(time,data.ydata,'ko','MarkerSize',1.5)
ax=gca;
ax.XLim = [min(time), max(time)];

end
