function Y=P3SIRoutput(tData,paramsNoK,Y0, outputSet)
paramsWithK=[paramsNoK(1) 1 paramsNoK(2) paramsNoK(3)];
ode_options = odeset('RelTol',1e-6);
[~,Y]=ode45(@SIR_rhs,tData,Y0,ode_options,paramsWithK);
if strcmpi(outputSet,'infected')
    Y=Y(:,2);
end
end

