function dy = SIR_rhs(t,y,params)
gamma = params(1);
k = params(2);
delta = params(3);
r = params(4);
dy = [delta*sum(y) - delta*y(1) - gamma*k*y(2)*y(1); %S
gamma*k*y(2)*y(1) - (r + delta)*y(2); %I
r*y(2) - delta*y(3)];%R
end

