%Project 5 main
clear; clc; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Problem 1
load q_data.txt;
q_data=sort(q_data);
M=15; k=8;
qm=2*rand(M,1);qm=sort(qm); 
% qm=linspace(0,2,M);
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
