function T = UninsulatedRodEquil(x,params)
%UNINSULATEDROD Summary of this function goes here
%   params:phi
Tambient=21.19;
a=.95;b=a;
L=70;

phi=params(1);
h=params(2);
k=params(3);

gamma=sqrt(2*(a+b)*h/(a*b*k));
kg=k*gamma;
c1=-phi/kg*((exp(gamma*L)*(h+kg))/(exp(-gamma*L)*(h-kg)+exp(gamma*L)*(h+kg)));
c2=phi/kg+c1;
T=c1*exp(-gamma*x)+c2*exp(gamma*x)+Tambient;
end
