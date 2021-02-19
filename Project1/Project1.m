%% Project 1 MA540
clear; close all
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
close all
k1=20.5; c1=1.5;
param1=[k1 c1];
h1=1e-16; % complex step
time=0:.01:6;

% Complex approx
kcomplex=complex(k1,h1);
paramsk=[kcomplex c1];
[t,ykcomplex]=ode15s(@(t,y) spring(t,y,paramsk), time, [2 ;-c1]); 
spring_k = imag(ykcomplex(:,1))/h1;

ccomplex=complex(c1,h1);
paramsc=[k1 ccomplex];
[t,yccomplex]=ode15s(@(t,y) spring(t,y,paramsc), time, [2 ;-ccomplex]); 
spring_c = imag(yccomplex(:,1))/h1;

% Finite Differences 
[t,yog]=ode15s(@(t,y) spring(t,y,param1), time, [2 ;-c1]); % original model for k
[t,yog2]=ode15s(@(t,y) spring(t,y,param1), time, [2 ;-1]); 
h2=1e-6;

% Finite k
parsfink=[k1+h2 c1];
[t,yk]=ode15s(@(t,y) spring(t,y,parsfink), time, [2 ;-c1]);
finitek=(yk(:,1)-yog(:,1))/h2;

% Finite c
parsfinc=[k1 c1+h2];
[t,yc]=ode15s(@(t,y) spring(t,y,parsfinc), time, [2 ;-1-h2]);
finitec=(yc(:,1)-yog2(:,1))/h2;

% Analytical solutions (given)
f0 = sqrt(4*k1 - c1^2);
f1 = sqrt(k1 - c1^2/4)*time;
f2 = exp(-c1*time/2);
f3 = (-2*time/f0).*sin(f1);
f4 = (c1*time/f0).*sin(f1) - time.*cos(f1);

kana=f2.*f3;
cana=f2.*f4;
ana=2*exp(-c1.*time./2).*cos(sqrt(k1-(c1^2)/4).*time);


figure('Renderer', 'painters', 'Position', [100 100 1050 550])
%k: currently all are the same
subplot(1,2,1)
hold on
plot(time,spring_k,LineSpec(1,'line'))
plot(time, finitek,LineSpec(2,'line'))
plot(time, kana,LineSpec(3,'line'))
hold off
legend('Complex','FD','Analytical','Location','Southeast')
title('Sensitivites of $k$ vs. Time','Interpreter','Latex')
xlabel('$Time$','Interpreter','Latex')
ylabel('Sensitivity, $\frac{\delta k}{\delta t}$', 'Interpreter','Latex')
axis([0 6 -0.5 0.5])

%c: currently all are same
% figure(2)
subplot(1,2,2)
hold on
plot(time,spring_c,LineSpec(1,'line'))
plot(time, finitec,LineSpec(2,'line'))
plot(time,cana,LineSpec(3,'line'))
hold off
legend('Complex','FD','Analytical','Location','Southeast')
title('Sensitivites of $c$ vs. Time','Interpreter','Latex')
xlabel('$Time$','Interpreter','Latex')
ylabel('Sensitivity, $\frac{\delta c}{\delta t}$', 'Interpreter','Latex')
axis([0 6 -0.5 0.5])

%Save Figure
saveas(gcf,'Figures/SpringSensitivities.png')

% figure(3) % Looking at ODE45 solution vs analytical solution
% hold on
% plot(time, ana)
% plot(time,yog(:,1))
% hold off

%% Problem 2

tf = 5;
dt = 0.01;
t_data = 0:dt:tf;
S0 = 900; R0 = 0; I0 = 100;
Y0 = [S0; I0; R0];
params = [0.2; 0.1; 0.15; 0.6];
ode_options = odeset('RelTol',1e-6);
[t,Y] = ode45(@SIR_rhs,t_data,Y0,ode_options,params);

% Finite differences to get Fisher matrix
deltah2=1e-6;

deltagamma=[0.2+deltah2;0.1;0.15;0.6];
[t,Ygam] = ode45(@SIR_rhs,t_data,Y0,ode_options,deltagamma);
pgam=(Ygam-Y)/deltah2; pgam=pgam(:,3);

deltak=[0.2; 0.1+deltah2; 0.15; 0.6];
[t,Yk] = ode45(@SIR_rhs,t_data,Y0,ode_options,deltak);
pk=(Yk-Y)/deltah2; pk=pk(:,3);

deltadelta=[0.2; 0.1; 0.15+deltah2; 0.6];
[t,Ydelta] = ode45(@SIR_rhs,t_data,Y0,ode_options,deltadelta);
pdelta=(Ydelta-Y)/deltah2; pdelta=pdelta(:,3);

deltar=[0.2; 0.1; 0.15; 0.6+deltah2];
[t,Yr] = ode45(@SIR_rhs,t_data,Y0,ode_options,deltar);
pr=(Yr-Y)/deltah2; pr=pr(:,3);

chi=[pgam pk pdelta pr];
fisher1=chi'*chi;
eig1=eig(fisher1);
[v,d]=eig(fisher1);

chi2=[pk pdelta pr];
fisher2=chi2'*chi2;
eig2=eig(fisher2);

%% Problem 3
%Set rod settings
rodPoints=10:4:70;
rodParams=[-18.4, .00191, 2.37];
%Set figure settings
hold on
for i=1:2
    if i==1
        figure('Renderer', 'painters', 'Position', [100 100 1050 550])
        %Get sensitivity matrices
        jacFinite=getJacobian(@(params)UninsulatedRodEquil(rodPoints,params),rodParams,'real');
        jacComplex=getJacobian(@(params)UninsulatedRodEquil(rodPoints,params),rodParams,'complex');
        for j=1:2
            if j==1
                subplot(1,2,1)
                for iParam=1:size(jacFinite,2)
                    semilogy(rodPoints,abs(jacComplex(:,iParam)),LineSpec(iParam,'line'))
                    hold on
                end
                ylabel('$\left|\frac{\partial T}{\partial \theta_i}\right|$','Interpreter','Latex')
                title({'Complex Sensitivities' ''},'Interpreter','Latex')
            elseif j==2
                subplot(1,2,2)
                jacAnalytic=RodSensitivities(rodParams,rodPoints');
                for iParam=1:size(jacComplex,2)
                    semilogy(rodPoints,abs(jacComplex(:,iParam)-jacAnalytic(:,iParam)),LineSpec(iParam,'line'))
                    hold on
                end
                title({'Complex - Analytic Difference' ''},'Interpreter','Latex')
                legend('$\phi=-18.4$','$h=.00191$','$k=2.37$','Interpreter','Latex')
            end
            xlabel('$x (cm)$','Interpreter','Latex')
        end
    elseif i==2
%         %Get sensitivity matrices
%         jacFinite=getJacobian(@(params)UninsulatedRodEquil(rodPoints,[params 2.37]),rodParams(1:2),'real');
%         jacComplex=getJacobian(@(params)UninsulatedRodEquil(rodPoints,[params 2.37]),rodParams(1:2),'complex');
%         for j=1:2
%             subplot(1,2,j)
%             hold on
%             ylabel('$\frac{\partial T_C}{\partial \theta_i}-\frac{\partial T_F}{\partial \theta_i}$','Interpreter','Latex')
%             if j==1
%                 for iParam=1:size(jacFinite,2)
%                     plot(rodPoints,jacComplex(:,iParam),LineSpec(iParam,'line'))
%                 end
%                 title({'Complex Sensitivities' ''},'Interpreter','Latex')
%             elseif j==2
%                 for iParam=1:size(jacFinite,2)
%                     plot(rodPoints,jacComplex(:,iParam)-jacFinite(:,iParam),LineSpec(iParam,'line'))
%                 end
%                 title({'Complex - Finite Difference' ''},'Interpreter','Latex')
%                 legend('$\phi=-18.4$','$h=.00191$','$k=2.37$','Interpreter','Latex')
%             end
%             xlabel('$x (cm)$','Interpreter','Latex')
%         end
    end
    
    %Plot Complex Sensitivity Results
    %Get fisher matrices
    scaledFisherFinite=jacFinite'*jacFinite;
    scaledFisherComplex=jacComplex'*jacComplex;

    %Compute ranks
    rankFiniteFisher=rank(scaledFisherFinite);
    rankComplexFisher=rank(scaledFisherComplex);
end

saveas(gcf,'Figures/RodSensitivities.png')
%% Problem 4
% Part A
alpha1=-389.4; alpha11=761.3; alpha111=61.5;
psi=@(P) alpha1*P.^2+alpha11*P.^4+alpha111*P.^6;
interval=-.8:.01:.8;
figure('Renderer', 'painters', 'Position', [100 100 800 600])
plot(interval, psi(interval),LineSpec(1,'line'))
title('Helmholtz Energy vs. Polarization')
xlabel('Polarization, P', 'Interpreter','Latex')
ylabel('Helmholtz Energy, $\psi$', 'Interpreter','Latex')
saveas(gcf,'Figures/HelmholtzEnergy.png')

% Part B
interval2=linspace(0,0.8,17);
deltalph1=interval2.^2;
deltalph2=interval2.^4;
deltalph3=interval2.^6;
chi4=[deltalph1' deltalph2' deltalph3'];
fisher4=chi4'*chi4;
eigenvalues4=eig(fisher4);

%Part C
%Generate Samples
sampNumberMorris=50;
alphaSample=[-389.4, 761.3, 61.5].*(rand(50,3)*.4+.8);     %Sample Unif(.8alpha,1.2alpha)
dMat=NaN(sampNumberMorris,3);
for i=1:sampNumberMorris
    dMat(i,:)=getJacobian(@(params)HelmholtzInt(params,[0 .8]), ...
        alphaSample(i,:),'finite','h',1/20);
end
muStarMorris=sum(abs(dMat),1)/sampNumberMorris;
muMorris=sum(dMat,1)/sampNumberMorris;
sigmaSMorris=1/(sampNumberMorris-1)*sum((dMat-muMorris).^2,1);
alphaJac=getJacobian(@(params)HelmholtzInt(params,[0 .8]),[-389.4, 761.3, 61.5]);
fprintf('\nProblem 4\n')
fprintf(['muStar=[%.3f %.3f %.3f]\nsigmaSquared=[%.3f %.3f %.3f]\n'...
    'alphaJac=[%.3f %.3f %.3f]\n'],muMorris,sigmaSMorris,alphaJac)


%Part D
sampNumberSobol=100000;
alphaRanges=[-389.4, 761.3, 61.5].*([.8;1.2].*ones(2,3));
[sobolBaseFull, sobolTotFull,sampleFull]=SatelliSobol(...
    @(params)HelmholtzInt(params,[0 0.8]),alphaRanges,sampNumberSobol);
    %param(3)=alpha111 found to be noninfluential and second order effects
    %minimal
[sobolBaseReduced, sobolTotReduced,sampleReduced]=SatelliSobol(...
    @(params)HelmholtzInt([params 61.5*ones(size(params,1),1)],[0 0.8]),alphaRanges(:,1:2),sampNumberSobol);

fprintf(['\nUnreduced Sobol: [%.3g, %.3g, %.3g]\n'...
         'Unreduced Total Sobol: [%.3g, %.3g, %.3g]\n'...
         'Reduced Sobol: [%.3g, %.3g]\n'...
         'Reduced Total Sobol: [%.3g, %.3g]\n'],...
         sobolBaseFull, sobolTotFull, sobolBaseReduced, sobolTotReduced)


[~,densityFull,xMeshFull,cdfFull]=kde(sampleFull,4000);
[~,densityReduced,xMeshReduced,cdfReduced]=kde(sampleFull,4000);
figure('Renderer', 'painters', 'Position', [100 100 800 600])
hold on
plot(xMeshFull,densityFull,LineSpec(1,'line'))
plot(xMeshReduced,densityReduced,LineSpec(2,'line'))
xlabel('','Interpreter','Latex')
ylabel('pdf','Interpreter','Latex')
title('Kernal Density Estimates of Heimholtz Energy','Interpreter','Latex')
legend('Unfixed $\alpha_{111}$','Fixed $\alpha_{111}$','Interpreter','Latex')
saveas(gcf,'Figures/HelmholtzKDE.png')


    




%% Function for Prob 1

function dy=spring(t,y,param1)
k=param1(1); c=param1(2);
dy(1)=y(2);
dy(2)=(-c*dy(1))-(k*y(1));
dy=[dy(1); dy(2)];
end

%% Function for Prob 2

function dy = SIR_rhs(t,y,params)
N = 1000;
gamma = params(1);
k = params(2);
delta = params(3);
r = params(4);
dy = [delta*N - delta*y(1) - gamma*k*y(2)*y(1); %S
gamma*k*y(2)*y(1) - (r + delta)*y(2); %I
r*y(2) - delta*y(3)];%R
end

%% Functions for Prob 3

%% Functions for Prob 4
function int=HelmholtzInt(params,yBounds)

int=[params(:,1)/3 params(:,2)/5 params(:,3)/7]*...
    [yBounds(2)^3-yBounds(1)^3; yBounds(2)^5-yBounds(1)^5; yBounds(2)^7-yBounds(1)^7];

end

function [sobolBase, sobolTot,sampC]=SatelliSobol(evalFcn,paramRange,sampleSize)
nPOIs=size(paramRange,2);
sampA=rand(sampleSize,size(paramRange,2)).*(paramRange(2,:)-paramRange(1,:))+paramRange(1,:);
sampB=rand(sampleSize,size(paramRange,2)).*(paramRange(2,:)-paramRange(1,:))+paramRange(1,:);
sampC=[sampA; sampB];
fA=evalFcn(sampA); fB=evalFcn(sampB); fC=evalFcn(sampC);
fAb=NaN(sampleSize,nPOIs);fBa=fAb;
for iParams=1:nPOIs
    sampAb=sampA; sampAb(:,iParams)=sampB(:,iParams);
    sampBa=sampB; sampBa(:,iParams)=sampA(:,iParams);
    fAb(:,iParams)=evalFcn(sampAb);
    fBa(:,iParams)=evalFcn(sampBa);
end
    fCexpected=mean(evalFcn(sampC));          %UNSURE IF THIS IS CORRECT
    sobolDen=1/(2*sampleSize)*sum(fC.^2,1)-fCexpected.^2;
    sobolBase=(1/sampleSize*sum(fA.*fBa-fA.*fB,1))/sobolDen;
    sobolTot=(1/(2*sampleSize)*sum(fA-fAb,1))/sobolDen;
end


