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
alphanom=[-389.4,761.3,61.5]; % nominal parameter estimates (reference)

psi=@(alpha,P) alpha(1)*P.^2+alpha(2)*P.^4+alpha(3)*P.^6; % function

res=udata-psi(alphaols,xdata); % residuals

design1=xdata.^2; design2=xdata.^4; design3=xdata.^6; % making design matrix
design=[design1 design2 design3]'; 

varhat=(1/(n-3))*(res'*res);
V = varhat.*inv(design*design');

%%%% DRAM
clear data model options

data.xdata = xdata;
data.ydata = udata;
tcov = V;
tmin = [alphaols(1); alphaols(2); alphaols(3)];

params = {
{'\alpha_1',tmin(1),-500}
{'\alpha_{11}',tmin(2),700}
{'\alpha_{111}',tmin(3),0}};

model.ssfun = @helmss;
model.sigma2 = varhat;
options.qcov = tcov;
options.nsimu = 10000;
options.updatesigma = 1;
N = 10000;
  
[results,chain,s2chain] = mcmcrun(model,data,params,options); % running DRAM

DRAM1vals = chain(:,1);
DRAM11vals = chain(:,2);
DRAM111vals=chain(:,3);

DRAMrange_1 = max(DRAM1vals) - min(DRAM1vals);
DRAMrange_11 = max(DRAM11vals) - min(DRAM11vals);
DRAMrange_111 = max(DRAM111vals) - min(DRAM111vals);
DRAM1_min = min(DRAM1vals)-DRAMrange_1/10;
DRAM1_max = max(DRAM1vals)+DRAMrange_1/10;
DRAM11_min = min(DRAM11vals)-DRAMrange_11/10;
DRAM11_max = max(DRAM11vals)+DRAMrange_11/10;
DRAM111_min = min(DRAM111vals)-DRAMrange_111/10;
DRAM111_max = max(DRAM111vals)+DRAMrange_111/10;
[bandwidth_DRAM1,density_DRAM1,DRAM1mesh,cdf_DRAM1]=kde(DRAM1vals);
[bandwidth_DRAM11,density_DRAM11,DRAM11mesh,cdf_DRAM11]=kde(DRAM11vals);
[bandwidth_DRAM111,density_DRAM111,DRAM111mesh,cdf_DRAM111]=kde(DRAM111vals);

cov(chain)
chainstats(chain,results)

figure('Renderer', 'painters', 'Position', [100 100 1050 700])
subplot(3,1,1)
plot(DRAM1vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_1$','Interpreter','Latex')

subplot(3,1,2);
plot(DRAM11vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_{11}$','Interpreter','Latex')

subplot(3,1,3);
plot(DRAM111vals,'-','linewidth',2)
set(gca,'Fontsize',[16]);
axis([0 N -inf inf])
xlabel('Chain Iteration')
ylabel('$\alpha_{111}$','Interpreter','Latex')
saveas(gcf,'Figures/EstimatesDRAMHelmholtz.png')

figure('Renderer', 'painters', 'Position', [100 100 1050 700])
subplot(1,3,1)
plot(DRAM1mesh,density_DRAM1,'k-','linewidth',3)
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')

subplot(1,3,2)
plot(DRAM11mesh,density_DRAM11,'k-','linewidth',3)
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_{11}$','Interpreter','Latex')

subplot(1,3,3)
plot(DRAM111mesh,density_DRAM111,'k-','linewidth',3)
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_{111}$','Interpreter','Latex')
saveas(gcf,'Figures/kdeDRAMHelmholtz.png')

figure('Renderer', 'painters', 'Position', [100 100 1050 700])
subplot(1,3,1)
scatter(DRAM1vals,DRAM11vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')
ylabel('$\alpha_{11}$','Interpreter','Latex')

subplot(1,3,2)
scatter(DRAM1vals,DRAM111vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_1$','Interpreter','Latex')
ylabel('$\alpha_{111}$','Interpreter','Latex')

subplot(1,3,3)
scatter(DRAM11vals,DRAM111vals)
box on
axis([-inf inf -inf inf])
set(gca,'Fontsize',[22]);
xlabel('$\alpha_{11}$','Interpreter','Latex')
ylabel('$\alpha_{111}$','Interpreter','Latex')
saveas(gcf,'Figures/CorrDRAMHelmholtz.png')

%%% RW Metropolis
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
%% Support Functions

function ss=sirSS(data,params,inits,useK)
if strcmpi(useK,'no k')
    params=[params(1) 1 params(2) params(3)];
end
ode_options = odeset('RelTol',1e-6);
[~,Y]=ode45(@SIR_rhs,data.xdata,inits,ode_options,params);
ss=sum((data.ydata-Y(:,2)).^2);
end
