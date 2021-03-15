%% Project 3
% Leah Rolf, Harley Hanes, James Savino
clc; clear all; close all

%% Problem 1

%% Problem 2

%% Problem 3

load Helmholtz.txt

xdata=Helmholtz(:,1); n=length(xdata); % allocating columns of data
udata=Helmholtz(:,2); 

alphaols=[-392.66,770.10,57.61]; % pulled from project 2 parameter estimates

psi=@(alpha,P) alpha(1)*P.^2+alpha(2)*P.^4+alpha(3)*P.^6; % function

res=udata-psi(alphaols,xdata); % residuals

design1=xdata.^2; design2=xdata.^4; design3=xdata.^6; % making design matrix
design=[design1 design2 design3]'; 

varhat=(1/(n-3))*(res'*res);
V = varhat.*inv(design*design');

% clear data model options % getting ready for DRAM
% 
% data.xdata = xdata';
% data.ydata = udata';
% tcov = V;
% tmin = [alphaols(1); alphaols(2); alphaols(3)];
% 
% params = {
% {'q1',tmin(1),-390}
% {'q2',tmin(2),770}
% {'q3',tmin(3),55};
% 
% model.ssfun = helmss(alphaols,data);
% model.sigma2 = varhat;
% options.qcov = tcov;
% options.nsimu = 10000;
% options.updatesigma = 1;
% N = 10000;
%   
% [results,chain,s2chain] = mcmcrun(model,data,params,options); % running DRAM

% RW Metropolis
N = 1e+5;
R = chol(V);
q_old = [alphaols(1);alphaols(2);alphaols(3)];
SS_old = res'*res;
n0 = 0.001;
sigma02 = varhat;
aval = 0.5*(n0 + 15);
bval = 0.5*(n0*sigma02 + SS_old);
sigma2 = 1/gamrnd(aval,1/bval);
accept = 0;
  
for i = 1:N
z = randn(3,1); 
q_new = q_old + R'*z;
alpha1 = q_new(1,1);
alpha11 = q_new(2,1);
alpha111=q_new(3,1);
res=udata-psi([alpha1 alpha11 alpha111],xdata);
SS_new = res'*res;
u_alpha = rand(1);
term = exp(-.5*(SS_new-SS_old)/sigma2);
alpha = min(1,term);
if u_alpha < alpha
  Q_MCMC(:,i) = [alpha1; alpha11;alpha111];
  q_old = q_new;
  SS_old = SS_new;
  accept = accept + 1;
else
  Q_MCMC(:,i) = q_old;
end
Sigma2(i) = sigma2;
bval = 0.5*(n0*sigma02 + SS_old);
sigma2 = 1/gamrnd(aval,1/bval);
end

acceptance=accept/N;
alpha1vals = Q_MCMC(1,:);
alpha11vals = Q_MCMC(2,:);
alpha111vals=Q_MCMC(3,:);

range_1 = max(alpha1vals) - min(alpha1vals);
range_11 = max(alpha11vals) - min(alpha11vals);
range_111 = max(alpha111vals) - min(alpha111vals);
alpha1_min = min(alpha1vals)-range_1/10;
alpha1_max = max(alpha1vals)+range_1/10;
alpha11_min = min(alpha11vals)-range_11/10;
alpha11_max = max(alpha11vals)+range_11/10;
alpha111_min = min(alpha111vals)-range_111/10;
alpha111_max = max(alpha111vals)+range_111/10;

[bandwidth_1,density_1,alpha1mesh,cdf_1]=kde(alpha1vals);
[bandwidth_11,density_11,alpha11mesh,cdf_11]=kde(alpha11vals);
[bandwidth_111,density_111,alpha111mesh,cdf_111]=kde(alpha111vals);

figure('Renderer', 'painters', 'Position', [100 100 1050 700])
subplot(3,1,1)
plot(alpha1vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_1$','Interpreter','Latex')

subplot(3,1,2)
plot(alpha11vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_{11}$','Interpreter','Latex')

subplot(3,1,3)
plot(alpha111vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_{111}$','Interpreter','Latex')
% saveas(gcf,'Figures/EstimatesRWHelmholtz.png')

figure(2)
plot(Sigma2)
set(gca,'Fontsize',[22]);
title('Measurement Error Variance \sigma^2')
% saveas(gcf,'Figures/ErrorRWHelmholtz.png')

figure(3)
subplot(1,2,1)
plot(alpha1mesh,density_1,'k-','linewidth',3)
axis([-inf inf 0 inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')

subplot(1,2,2)
plot(alpha11mesh,density_11,'k-','linewidth',3)
axis([-inf inf 0 inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_{11}$','Interpreter','Latex')
% saveas(gcf,'Figures/kdeRWHelmholtz.png')

figure('Renderer', 'painters', 'Position', [100 100 1050 700])
subplot(1,3,1)
scatter(alpha1vals,alpha11vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')
ylabel('$\alpha_{11}$','Interpreter','Latex')

subplot(1,3,2)
scatter(alpha1vals,alpha111vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')
ylabel('$\alpha_{111}$','Interpreter','Latex')

subplot(1,3,3)
scatter(alpha11vals,alpha111vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_{11}$','Interpreter','Latex')
ylabel('$\alpha_{111}$','Interpreter','Latex')
% saveas(gcf,'Figures/CorrRWHelmholtz.png')