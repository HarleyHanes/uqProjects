%%
%                  Example16_15.m 
%
%
%  Compute the mean and standard deviation of the steady state spring solution
%
% $$ y(\omega_F,Q) = \frac{1}{\sqrt((k-m \omega_F^2)^2 + (c\omega_F)^2} $$
%
%  with coefficients Q = [m,c,k].  Here we compare three methods:
%   * Monte Carlo sampling as discussed in Section 13.3 with mean and variance at omega_f given by (13.23)
%   * Deterministic solutions y(omega_F,\bar{q})
%   * Discrete projection
%
%
  clear all

%
%

  K = 4;          % Number of tensored Legendre basis functions Psi(q)
  N_MC = 1e5;      % Number of Monte Carlo samples
  M = 10;          % Number of Gauss-Legendre quadrature points used to approximate discrete projection

%%
% Input the parameter means and standard deviations.
%


  mu_k = 8.5;
  sigma_k = 0.001;

%%
% Compute N_MC samples for each of the parameters and construct the grid for omega_F.
%

  q = zeros(N_MC,1);
  q = (mu_k - sigma_k) + 2*sigma_k*rand(N_MC,1);
 
  d_t = 0.01;                   % d_omega = 0.001; Example16_10.m 
  pt1 = 0;
  pt2 = 2.15;
  time = [pt1:d_t:pt2];
  N_time = length(time);

%%
% Input the M=10 Gauss-Legendre quadrature points and weights from a table.  Recall that one
% must multiply the weights by 0.5 to account for the density \rho(q) = 1/2. 
%

  zvals = [-0.9739 -0.86506 -0.6794 -0.4334 -0.14887 0.14887 0.4334 0.6794 0.86506 0.9739];
  weights = 0.5*[0.06667 0.14945 0.21908 0.2693 0.2955 0.2955 0.2693 0.21908 0.14945 0.06667];
  
  sigmak = 2*sigma_k/sqrt(12);
  kvals = mu_k*ones(size(zvals)) + sqrt(3)*sigmak*zvals;

%%
% Constuct the Legendre polynomials P_0, P_1 and P_2 evaluated at the quadrature points zvals. 
% Also construct the coefficients h0, h1 and h2 used to construct gamma_k.
%

  poly0 = ones(1,M);                              % P_0(q) = 1
  poly1 = zvals;                                  % P_1(q) = q
  poly2 = (3/2)*zvals.^2 - (1/2)*ones(1,M);       % P_2(q) = (3/2)q^2 - 1/2
  h0 = 1;
  h1 = 1/3;
  h2 = 1/5;

  y_k0 = zeros(1,N_time);
  y_k1 = zeros(1,N_time);
  y_k2 = zeros(1,N_time);
  y_k3 = zeros(1,N_time);

%%
%  Compute the deterministic, Monte Carlo and discrete projection values at each of values omega_F.
%

  for k = 1:N_time
    deterministic(k) = 3*cos(sqrt(mu_k)*time(k));     % Deterministic solution
    MC = 3*cos(sqrt(q(:,1))*time(k));          % Monte Carlo samples
    MC_mean(k) = (1/N_MC)*sum(MC);                                               % MC mean
    tmp = 0;
    for j=1:N_MC
       tmp = tmp + (MC(j) - MC_mean(k))^2;
    end
    MC_sigma(k) = sqrt((1/(N_MC-1))*tmp);

    for j1 = 1:M
      for j2 = 1:M
        for j3 = 1:M
          yval = 3*cos(sqrt(kvals(j1))*time(k));
          y_k0(k) = y_k0(k) + (1/(h0*h0*h0))*yval*poly0(j1)*poly0(j2)*poly0(j3)*weights(j1)*weights(j2)*weights(j3);
          y_k1(k) = y_k1(k) + (1/(h1*h0*h0))*yval*poly1(j1)*poly0(j2)*poly0(j3)*weights(j1)*weights(j2)*weights(j3);
          y_k2(k) = y_k2(k) + (1/(h0*h1*h0))*yval*poly0(j1)*poly1(j2)*poly0(j3)*weights(j1)*weights(j2)*weights(j3);
          y_k3(k) = y_k3(k) + (1/(h0*h0*h1))*yval*poly0(j1)*poly0(j2)*poly1(j3)*weights(j1)*weights(j2)*weights(j3);
        end
      end 
    end
  end

%%
% Compute the mean and variance based on the generalized Fourier coefficients y_k.
%

  Y_k = [y_k0; y_k1; y_k2; y_k3];
  gamma = [h0*h0*h0 h1*h0*h0 h0*h1*h0 h0*h0*h1];

  for n = 1:N_time
    var_sum = 0;
    for j = 2:K
      var_sum = var_sum + (Y_k(j,n)^2)*gamma(j);
    end
    var(n) = var_sum;
  end

  DP_mean = Y_k(1,:);           % Discrete projection mean
  DP_sd = sqrt(var);            % Discrete projection standard deviation
%%
%  Plot the response mean and standard deviation as a function of the drive frequency omega_F.
%

figure('Renderer', 'painters', 'Position', [100 100 1920 1080])
%  plot(omega,DP_mean,'b-',omega,MC_mean,'--r',omega,deterministic,'-.k','linewidth',3)
  plot(time,DP_mean,'k-',time,MC_mean,'--m','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('time')
  ylabel('Response Mean')
%  legend('Discrete Projection','Monte Carlo','Deterministic Solution','Location','NorthEast')
  legend('Discrete Projection','Monte Carlo','Location','North')
  saveas(gcf,'Figures/P2_Means')

  
figure('Renderer', 'painters', 'Position', [100 100 1920 1080])
  plot(time,DP_sd,'-k',time,MC_sigma,'--m','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('time')
  ylabel('Response Standard Deviation')  
  legend('Discrete Projection','Monte Carlo','Location','NorthWest')
  saveas(gcf,'Figures/P2_STD')
  
  
figure('Renderer', 'painters', 'Position', [100 100 1920 1080])
%  plot(omega,DP_mean,'b-',omega,MC_mean,'--r',omega,deterministic,'-.k','linewidth',3)
  plot(time,DP_mean-MC_mean,'k-','linewidth',3)
  %axis([pt1 pt2 -3 3])
  set(gca,'Fontsize',[22]);
  xlabel('time')
  ylabel('Response Mean difference')
%  legend('Discrete Projection','Monte Carlo','Deterministic Solution','Location','NorthEast')
  %legend('Discrete Projection','Monte Carlo','Location','NorthEast')
  saveas(gcf,'Figures/P2_meanDiff')
   
figure('Renderer', 'painters', 'Position', [100 100 1920 1080])
%  plot(omega,DP_mean,'b-',omega,MC_mean,'--r',omega,deterministic,'-.k','linewidth',3)
  plot(time,DP_sd-MC_sigma,'k-','linewidth',3)
  %axis([pt1 pt2 -3 3])
  set(gca,'Fontsize',[22]);
  xlabel('time')
  ylabel('Response std difference')
%  legend('Discrete Projection','Monte Carlo','Deterministic Solution','Location','NorthEast')
  %legend('Discrete Projection','Monte Carlo','Location','NorthEast')
  saveas(gcf,'Figures/P2_stdDiff')

