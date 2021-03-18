%% Project 3
% Leah Rolf, Harley Hanes, James Savino
clc; clear all; close all
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1

%% Problem 2
    problem='4 parameter';
    %Initialize problem
        %Load SIR data
        D=load('SIR.txt');
        data.xdata=D(:,1);
        data.ydata=D(:,2);
        clear('D')
        %Set initial conditions for SIR
        inits=[1000-127.1233; 127.1233; 0];
        
        %Seperate part A and B
        if strcmpi(problem,'3 parameter')
            %State parameter names and ranges
            params={
                    {'gamma',0.00921,0,1}
                    {'delta', .195,0,1}
                    {'r', .787,0,1}
                };
            %Define error function
            model.ssfun=@(params,data)sirSS(data,params,inits,'no k');
            %Define parameter variances
            options.qcov=[3.996604557707984e-08,7.294589132883216e-08,2.080488101816833e-07;7.294589132883216e-08,9.379222001817734e-05,9.841326827550314e-05;2.080488101816833e-07,9.841326827550314e-05,1.720799701083470e-04];
        else
            %State parameter names and ranges
            params={
                    {'gamma',0.1294,0,1}
                    {'k',.1163,0,1}
                    {'delta', .1948,0,1}
                    {'r', .7834,0,1}
                };
            %Define Error Function
            model.ssfun=@(params,data)sirSS(data,params,inits,'with k');
        end
        %Set observation variance
        model.sigma2=3.890110168320449e+02;
        %Set number of simulations
        options.nsimu=10000;
        options.updatesigma=1;
        
    %Run DRAM
        [results,chain,s2chain]=mcmcrun(model,data,params,options);
    %Print Chain results
        cov(chain)
        chainstats(chain,results)
        if strcmpi(problem,'3 parameter')
            fprintf('\nsigma^2 mean=%.4e\n',mean(s2chain))
            fprintf('V_{DRAM}-V_{OLS}')
            results.qcov-options.qcov
            fprintf('(V_{DRAM}-V_{OLS})./V_{DRAM}')
            (results.qcov-options.qcov)./results.qcov
        end
    %Get distribution
    if strcmpi(problem,'3 parameter')
      [bandwidth_g,density_g,gMesh,cdf_q]=kde(chain(:,1));
      [bandwidth_d,density_d,dMesh,cdf_d]=kde(chain(:,2));
      [bandwidth_r,density_r,rMesh,cdf_r]=kde(chain(:,3));
    else 
      [bandwidth_g,density_g,gMesh,cdf_q]=kde(chain(:,1));
      [bandwidth_k,density_k,kMesh,cdf_k]=kde(chain(:,2));
      [bandwidth_d,density_d,dMesh,cdf_d]=kde(chain(:,3));
      [bandwidth_r,density_r,rMesh,cdf_r]=kde(chain(:,4));
    end
        
    %Plot Results
        
%         figure; clf
%         mcmcplot(chain,[],results,'chainpanel');
        %Pairwise plots
        figure('Renderer', 'painters', 'Position', [100 100 650 550]); clf
        mcmcplot(chain,[],results,'pairs');
        saveas(gcf,sprintf('Figures/%s_pairwise',problem))

        %Parameter Chains
        if strcmpi(problem,'3 parameter')
        for i=1:3
            figure('Renderer', 'painters', 'Position', [100 100 650 450]); clf
            plot(chain(:,i),'-','LineWidth',1)
            %set(gca,'Fontsize',[22]);
            %axis([0 options.n -19.3 -17.5])
            xlabel('Chain Iteration')
            switch i
                case 1
                    ylabel(sprintf('$\\gamma\\ (\\bar{\\gamma}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_gChain',problem))
                case 2
                    ylabel(sprintf('$\\delta\\ (\\bar{\\delta}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_dChain',problem))
                case 3
                    ylabel(sprintf('$r\\ (\\bar{r}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_rChain',problem))
            end
        end
        else
        for i=1:4
            figure('Renderer', 'painters', 'Position', [100 100 650 450]); clf
            plot(chain(:,i),'-','LineWidth',1)
            %set(gca,'Fontsize',[22]);
            %axis([0 options.n -19.3 -17.5])
            xlabel('Chain Iteration')
            switch i
                case 1
                    ylabel(sprintf('$\\gamma\\ (\\bar{\\gamma}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_gChain',problem))
                case 2
                    ylabel(sprintf('$k\\ (\\bar{k}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_kChain',problem))
                case 3
                    ylabel(sprintf('$\\delta\\ (\\bar{\\delta}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_dChain',problem))
                case 4
                    ylabel(sprintf('$r\\ (\\bar{r}=%.4g)$',mean(chain(:,i))),'Interpreter','Latex')
                    saveas(gcf,sprintf('Figures/%s_rChain',problem))
            end
        end
        end
        %sigma squared chain
        figure('Renderer', 'painters', 'Position', [100 100 700 500])
        plot(s2chain,'-','LineWidth',1)
        xlabel('Chain Iteration')
        ylabel(sprintf('$\\sigma^2(\\bar{\\sigma^2}=%.4g)$',mean(s2chain)),'Interpreter','Latex')
        fprintf('\nsigma mean=%.4e\n',mean(s2chain))
        saveas(gcf,sprintf('Figures/%s_s2chain',problem))
        
        %Parameter distributions
        %gamma
        figure('Renderer', 'painters', 'Position', [100 100 650 450])
        plot(gMesh,density_g)
        xlabel('$\gamma$','Interpreter','Latex')
        saveas(gcf,sprintf('Figures/%s_Marginal_g',problem))
        
        %delta
        figure('Renderer', 'painters', 'Position', [100 100 650 450])
        plot(dMesh,density_d)
        xlabel('$\delta$','Interpreter','Latex')
        saveas(gcf,sprintf('Figures/%s_Marginal_d',problem))
        
        %k
        figure('Renderer', 'painters', 'Position', [100 100 650 450])
        plot(rMesh,density_r)
        xlabel('$k$','Interpreter','Latex')
        saveas(gcf,sprintf('Figures/%s_Marginal_k',problem))
        
        %h
        if ~strcmpi('3 parameter',problem)
        figure('Renderer', 'painters', 'Position', [100 100 650 450])
        plot(kMesh,density_k)
        xlabel('$h$','Interpreter','Latex')
        saveas(gcf,sprintf('Figures/%s_Marginal_h',problem))
        end

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

%% Support Functions

function ss=sirSS(data,params,inits,useK)
if strcmpi(useK,'no k')
    params=[params(1) 1 params(2) params(3)];
end
ode_options = odeset('RelTol',1e-6);
[~,Y]=ode45(@SIR_rhs,data.xdata,inits,ode_options,params);
ss=sum((data.ydata-Y(:,2)).^2);
end