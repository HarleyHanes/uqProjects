function [paramsOptimal,residualsOptimal,varEst,covEst] = LSQnonlin(evalFcn,data,problem)
%LSQnonlin Calculated optimal parameter values using nonlinear linear least squares
%Inputs
    %evalFcn(data.x,params):function oupting nx1 response vector y from nx1
    %   input vector x and px1 vector params
    %data: struct with .x and .y nx1 elements for input and output data
    %paramInit: initial parameter estimate for fminsearch
%Outputs
    %paramsOptimal: px1 Optimal parameter values
    %residualsOptimal:nx1 of residuals at optimal parameter value
    %varEst: scalar for estimate for sigma^2
    %covEst: pxp matrix for 
    
    %Check that data.x and data.y are column vectors of same length
    data=CheckDataSizes(data);
    %Make Cost function
    problem.objective=@(params)sum(sum((evalFcn(data.x,params)-data.y).^2));
    %Find Params with fminsearc
    paramsOptimal=fmincon(problem);
    %Get Residuals
    residualsOptimal=data.y-evalFcn(data.x,paramsOptimal);
    %Get variance estimate
    varEst=1/(length(data.x)-length(paramsOptimal))*(residualsOptimal'*residualsOptimal);
    %Get Covariance Estimate
        %Get Jacobian
        X=getJacobian(@(params)evalFcn(data.x,params),paramsOptimal);
    covEst=varEst*(X'*X)^(-1);
end

