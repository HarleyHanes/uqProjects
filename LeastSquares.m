function [paramsOptimal,rssOptimal,residualsOptimal] = LeastSquares(evalFcn,data,nParams,type)
%LeastSquares Calculated optimal parameter values using linear or nonlinear
%least squares
%Inputs
    %evalFcn(data.x,params):function oupting nxp matrix A where
    %Aij=params(j)*basefunction(j)*data.x(i)
    %data: struct with .x and .y elements for input and output data
    %       -needs to be oriented where each row is data entry
    %nParams: # of model params
    %type: 'linear' or 'nonlinear' for Least Squares method
%Outputs
    %paramsOptimal: Optimal parameter values
    %rssOptimal: Residual Sum of Squares at optimal parameter value
    %residualsOptimal:Vector of residuals at optimal parameter value
    baseParams=ones(1,nParams);
if strcmpi(type,'linear')
    %Make information matrix
    X=evalFcn(data.x,baseParams);
    %Caculate optimal parameters with pseudo-inverse
    paramsOptimal=(X'*X)^(-1)*X'*data.y;
    %Get Residuals
    residualsOptimal=data.y-sum(evalFcn(data.x,paramsOptimal),2);
elseif strcmpi(type,'nonlinear')
    %Make Cost function
    costFcn=@(params)sum(sum((evalFcn(data.x,ones(1,baseParams))-Data.y).^2));
    %Find Params with fminsearc
    paramsOptimal=fminsearch(costFcn,coeff0);
    %Get Residuals
    residualsOptimal=data.y-evalFcn(data.x,paramsOptimal);
else
    error('unrecognized Least Squares Type')
end
rssOptimal=sum(residualsOptimal.^2);
end

