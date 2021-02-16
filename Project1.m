%% Leah Rolf Project 1 MA540
clear all
close all

%% Problem 1
k1=20.5; c1=1.5;
param1=[k1 c1];
h1=1e-16;
time=1:.1:6;

% Making each variable complex
kcomplex=complex(k1,h1);
paramsk=[kcomplex c1];
[t,y]=ode45(@(t,y) spring(t,y,paramsk), time, [2 ;-c1]); % want y(:,2)
spring_k = imag(y(:,2))/h1;

ccomplex=complex(c1,h1);
paramsc=[k1 ccomplex];
[t,y]=ode45(@(t,y) spring(t,y,paramsc), time, [2 ;-c1]); % want y(:,2)
spring_c = imag(y(:,2))/h1;

% Analytical solutions (given)
kana=exp(-c1.*time./2).*((-2.*time)/sqrt(4*k1-(c1)^2)).*sin(time.*sqrt(k1-(c1^2)/4));
cana=exp(-c1.*time./2).*(c1.*time)/sqrt(4*k1-(c1)^2).*sin(time.*sqrt(k1-(c1^2)/4))-time.*cos(sqrt(k1-(c1^2)/4).*time);

% Finite Differences 
[t,y]=ode45(@(t,y) spring(t,y,param1), time, [2 ;-c1]); % original model
h=1e-6;

% Finite k
parsfink=[k1+h c1];
[t,yk]=ode45(@(t,y) spring(t,y,parsfink), time, [2 ;-c1]);
finitek=(yk(:,2)-y(:,2))/h;

% Finite h
parsfinc=[k1 c1+h];
[t,yh]=ode45(@(t,y) spring(t,y,parsfinc), time, [2 ;-c1]);
finitec=(yh(:,2)-y(:,2))/h;

figure(1) %currently complex and finite are the same
hold on
plot(time,spring_k)
plot(time, kana)
plot(time, finitek)
hold off

figure(2) %currently complex and finite are the same
hold on
plot(time,spring_c)
plot(time,cana)
plot(time, finitec)
hold off

%% Problem 2

tf = 5;
dt = 0.01;
t_data = 0:dt:tf;
S0 = 900; R0 = 0; I0 = 100;
Y0 = [S0; I0; R0];
params = [0.2; 0.1; 0.15; 0.6];
ode_options = odeset('RelTol',1e-6);
[t,Y] = ode45(@SIR_rhs,t_data,Y0,ode_options,params);

% Finite differences
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

%% Problem 4
Part A
alpha1=-389.4; alpha11=761.3; alpha111=61.5;
psi=@(P) alpha1*P.^2+alpha11*P.^4+alpha111*P.^6;
interval=-.8:.01:.8;
figure(20)
plot(interval, psi(interval), 'LineWidth', 3)
title('Helmholtz Energy vs. Polarization', 'FontSize',16)
xlabel('Polarization, P', 'Interpreter','Latex', 'FontSize',12)
ylabel('Helmholtz Energy, $\psi$', 'Interpreter','Latex', 'FontSize',12)

% Part B
interval2=linspace(0,0.8,17);
deltalph1=interval2.^2;
deltalph2=interval2.^4;
deltalph3=interval2.^6;
chi4=[deltalph1' deltalph2' deltalph3'];
fisher4=chi4'*chi4;
eigenvalues4=eig(fisher4);


%% Function for Prob 1

function dy=spring(t,y,param1)
dy(1)=y(2);
dy(2)=-param1(2)*y(2)-param1(1)*y(1);
dy=[dy(1); dy(2)];
end

%% Function for Prob 2

function dy = SIR_rhs(t,y,params)
N = 1000;
gamma = params(1);
k = params(2);
delta = params(3);
r = params(4);
dy = [delta*N - delta*y(1) - gamma*k*y(2)*y(1);
gamma*k*y(2)*y(1) - (r + delta)*y(2);
r*y(2) - delta*y(3)];
end
