%Project 5 main
clear; clc; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Problem 3
fEval=@(q)(6*q.^2+3).*sin(6*q-4);
xs=[3.8055214e-01
      1.7092556e-01
   6.2833659e-01
   7.1959107e-01
   1.9375821e+00
   1.6575153e+00
   9.4462033e-01
   1.4111168e+00
   5.2149567e-01
   3.1004611e-02
   1.2837187e+00
   8.6819959e-01
   1.8276181e+00
   1.1183761e+00
   1.4867790e+00];
ys=fEval(xs);
xpIn=linspace(0,2,100)';
xpOut=linspace(-.5,2.5,100)';
% 
% plot(qData,'*')
% hold on
% x=linspace(0,2,100);
% plot(x,feval(x))

%gpModel=fitrgp(qData,fEval(qData),'Basis','linear',...
%      'FitMethod','exact','PredictMethod','exact');
gprMdl = fitrgp(xs,ys,'KernelFunction','squaredexponential','KernelParameters',[3.5, 10],'Sigma',eps);

%ypred = resubPredict(gpModel)
[predIn,~,yintIn] = predict(gprMdl,xpIn);
[predOut,~,yintOut] = predict(gprMdl,xpOut);


figure('Renderer', 'painters', 'Position', [100 100 600 500])
  f = [yintIn(:,2); flipud(yintIn(:,1))];
  h(1) = fill([xpIn; flipud(xpIn)], f, [7 7 7]/8);
  set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  hold on
  h(2) = plot(xs,ys,'ro','linewidth',5,'DisplayName','Data');
  h(3) = plot(xpIn,predIn,'b-','linewidth',1,'DisplayName','Predictive Mean');
  hold off
  legend('Location','NorthWest')
  set(gca,'Fontsize',[22]);
  xlabel('Parameter q','Interpreter','Latex')
  ylabel('Response','Interpreter','Latex')
saveas(gcf,'Figures/Problem3_InPredict.png')
  
figure('Renderer', 'painters', 'Position', [100 100 600 500])
  f = [yintOut(:,2); flipud(yintOut(:,1))];
  h(1) = fill([xpOut; flipud(xpOut)], f, [7 7 7]/8);
  set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  hold on
  h(2) = plot(xs,ys,'ro','linewidth',5,'DisplayName','Data');
  h(3) = plot(xpOut,predOut,'b-','linewidth',1,'DisplayName','Predictive Mean');
  hold off
  legend('Location','NorthWest')
  set(gca,'Fontsize',[22]);
  xlabel('Parameter q','Interpreter','Latex')
  ylabel('Response','Interpreter','Latex')
  axis([-.5 2.5 -30 30])
saveas(gcf,'Figures/Problem3_OutPredict.png')
