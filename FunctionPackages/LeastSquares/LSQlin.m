function [paramsOptimal,residualsOptimal,varEst,covEst] = LSQlin(evalFcn,data)
%LSQlin Calculated optimal parameter values using linear least squares
%Inputs
    %evalFcn(data.x,params):function oupting nxp matrix A where
        %Aij=basefunction_j(data.x_i)
    %data: struct with .x and .y nx1 elements for input and output data
%Outputs
    %paramsOptimal: px1 Optimal parameter values
    %residualsOptimal:nx1 of residuals at optimal parameter value
    %varEst: scalar for estimate for sigma^2
    %covEst: pxp matrix for 
    
    %Check that data.x and data.y are column vectors of same length
    data=CheckDataSizes(data);
    %Make design matrix
    X=evalFcn(data.x);      
    %Caculate optimal parameters with pseudo-inverse
    paramsOptimal=(X'*X)\X'*data.y;
    %Get Residuals
    residualsOptimal=data.y-evalFcn(data.x)*paramsOptimal;
    %Calculate variance estimate
    varEst=1/(length(data.x)-length(paramsOptimal))*(residualsOptimal'*residualsOptimal);
    %Calculate Covariance Estimate
    covEst=varEst*(X'*X)^(-1);
end

