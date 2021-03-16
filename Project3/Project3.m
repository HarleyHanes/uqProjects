%% Project 3
% Leah Rolf, Harley Hanes, James Savino
clc; clear all; close all

%% Problem 1

%% Problem 2
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
            {'gamma',0.02,0,1}
            {'delta', .15,0,1}
            {'r', .6,0,1}
        };
        %Define error func
        model.ssfun=@(params,data)sirSS(data,params,inits,'no k');
        %Define varainces- taken from project 2 problem 3
        model.sigma2=3.890110168320449e+02;
        options.qcov=[3.996604557707984e-08,7.294589132883216e-08,2.080488101816833e-07;7.294589132883216e-08,9.379222001817734e-05,9.841326827550314e-05;2.080488101816833e-07,9.841326827550314e-05,1.720799701083470e-04];
        %Set number of simulations
        options.nsimu=10000;
        options.updatesigma=1;
        
    %Run DRAM
        [results,chain,s2chain]=mcmcrun(model,data,params,options);
        
    %Plot Results
        figure(1); clf
        mcmcplot(chain,[],results,'chainpanel');

        figure(2); clf
        mcmcplot(chain,[],results,'pairs');

        cov(chain)
        chainstats(chain,results)
        
        for i=1:3
            figure(2+i); clf
            plot(chain(:,i),'-')
            set(gca,'Fontsize',[22]);
            %axis([0 options.n -19.3 -17.5])
            xlabel('Chain Iteration')
            switch i
                case 1
                    ylabel('Parameter gamma')
                case 2
                    ylabel('Parameter delta')
                case 3
                    ylabel('Parameter r')
            end
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
    paramsWithK=[params(1) 1 params(2) params(3)];
end
ode_options = odeset('RelTol',1e-6);
[~,Y]=ode45(@SIR_rhs,data.xdata,inits,ode_options,paramsWithK);
ss=sum((data.ydata-Y(:,2)).^2);
end
