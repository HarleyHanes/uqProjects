function jac = getJacobian(evalFcn,xBase,varargin)
% Takes a set of parameters for the function bHandle
%  as well as any arguments you'd like to pass to bHandle
%
% POIs - baseline parameters, p
% evalFcn - q(params)
%
% Returns the normalized jacobian evaluated at p:
%   (pj / qi) (partial qi/ partial pj) = S

%%
if sum(strcmpi(varargin,'finite'))
    approxMethod='finite';
elseif sum(strcmpi(varargin,'forward'))
    approxMethod='forward';
else 
    approxMethod='complex';
end
if isempty(varargin)||(sum(strcmpi(varargin,'h'))==0) 
    delta_x0 = 1e-8;
elseif sum(strcmpi(varargin,'h'))==1
    delta_x0=varargin{find(strcmpi(varargin,'h'))+1};      %Assign whatever follows 'h' to delta_x0
end
nPOIs=length(xBase);
% calculate baseline point
baseQuants = evalFcn(xBase);

nQuants = length(baseQuants);
%raw=1;

%%
jac = NaN(nQuants, nPOIs);

%factor = .001; % get .1% of initial parameters, HARD

% fudge number if delta is 0
%delta_x0(delta_x0==0) = factor;
y0 = baseQuants;
for i = 1:nPOIs
    xPert=xBase;
    if strcmpi(approxMethod,'finite')
        % basic centered difference approximation of Jacobian
        xPert(i) = xBase(i) - delta_x0;
        yLo = evalFcn(xPert);
        
        xPert(i) = xBase(i) + delta_x0;
        yHi = evalFcn(xPert);
    elseif strcmpi(approxMethod,'complex')
        xPert(i)=xBase(i)+delta_x0*1i;
        yPert=evalFcn(xPert);
    elseif strcmpi(approxMethod,'forward')
        xPert(i)=xBase(i)+delta_x0;
        yPert=evalFcn(xPert);
    end
    for j = 1:nQuants
        % Calculate partial
        if strcmpi(approxMethod,'finite')
            jac(j,i) = (yHi(j) - yLo(j)) / (2 * delta_x0);
        elseif strcmpi(approxMethod,'complex')
            jac(j,i)=imag(yPert(j))/delta_x0;
        elseif strcmpi(approxMethod,'forward')
            jac(j,i)=(yPert-baseQuants)/delta_x0;
        end
       
%         if ~raw
%             jac(j,i) = jac(j,i) * init_x0(i) * sign(y0q) / sqrt(eps + y0q^2);
%         end
    end
end

% for iPOI=1:nPOIs
%         xMinus=x0(iPOI)-xDelta;
%         xPlus=x0(iPOI)+xDelta;
%         yMinus=evalFcn(xMinus);
%         yPlus=evalFcn(xPlus);
%     for iQuant=1:nQaunts
%         y0=
%         jac(iQuant,iPOI)=(yMinus(iQuant)-2*y0(iQuant,iPOI)+yPlus(iQuant))/(xDelta^2);
%     end
%         
% end


end