function jac=RodSensitivities(params,x)
    %Apply parameters
        Tamb=21.19;
        a=.95;b=a;
        L=70;
    %Get Symbolic Sensitivities
    syms h phi k
    Gamma=sqrt(2*(a+b)*h/(a*b*k));
    c1=-phi/(k*Gamma)*...
        (exp(Gamma*L)*(h+k*Gamma)/...
        (exp(-Gamma*L)*(h-k*Gamma)+exp(Gamma*L)*(h+k*Gamma)));
    c2=phi/(k*Gamma)+c1;
    y=c1*exp(-Gamma.*x)+c2*exp(Gamma.*x)+Tamb;

    kSensitivity=diff(y,k);
    phiSensitivity=diff(y,phi);
    hSensitivity=diff(y,h);
    
    jac=[phiSensitivity hSensitivity kSensitivity];

    
    jac=subs(jac,phi,params(1));
    jac=subs(jac,h,params(2));
    jac=subs(jac,k,params(3));
    
end