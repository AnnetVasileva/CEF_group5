%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is our Code based on...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Matlab
clc;
clear all;
close all;
%
%

% Amount of data points
long	= 20000;
init	= 500;
slong	= init+long;
% Define sequence in order to compute capital change in the end
T		= init+1:slong-1;
T1		= init+2:slong;
% Tolerance criterium for the loop
tol	= 1e-6;
% Condition for the main loop
crit	= 1;
% Updating parameter
gam	= .25;
% Risk-Aversion (Arrow-Pratt)
sigma	= 1;
% Depreciation rate
delta	= 0.1;
% Discount Factor
beta	= 0.95;
% Capital share from output
alpha	= 0.3;
ab		= 0;
% Persistence of shock
rho	= 0.8;
% Standard deviation for noise
se		= 0.125;

% In the following we compute the Steady State for Capital and Consumption
% Steady State of capital
kss	= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
% Steady State of consumption
css	= kss^alpha-delta*kss;

% Set seed for reproducibility
randn('state',1);
% Generate the shock
e		= se*randn(slong,1);
% Allocate memory for the simulated shock
a		= zeros(slong,1);
% We start with a shock of 0
a(1)	= ab+e(1);
% Generate the productivity shock
for i	= 2:slong;
    %a(i)=0.85*a(i-1)+0.15*0+e
   a(i)=rho*a(i-1)+(1-rho)*ab+e(i);
end

% Initialize the Algorithm Parameters
param	= [ab alpha beta delta rho se sigma long init];
% We want to start with value in steady state
% Our economy is in the balanced growth path (BGP)
b0		= peaoginit_with_explanation(e,param);

%
% Main Loop
%
iter		= 1;
while crit>tol;
   %
   % Allocate memory for the simulated serie
   % Capital stock
   k		= zeros(slong+1,1);
   % Lambda
   lb		= zeros(slong,1);
   % Lagrange Multiplier
   mu		= zeros(slong,1);
   % Generate Memory Matrice for State Variable (k and a(i))
   X		= zeros(slong,length(b0));
   % We set initial value of capital to steady state
   k(1)	= kss;
   for i	= 1:slong;
       % X is the Matrice with our State Variable
       % lb is our function with which we make the guess
       % Reason: agent's conditional expectations are time functions of models' state
       % variable k and a(i)
       % To account for the numerical instability we leave out the last part as
       % Marcet Hannes[1999]log(k(i))*a(i)
       
      X(i,:)= [1 log(k(i)) a(i) log(k(i))*log(k(i)) a(i)*a(i)];
      % lambda
      lb(i)	= exp(X(i,:)*b0);
      % Compute investment
      iv		= exp(a(i))*k(i)^alpha-lb(i)^(-1/sigma);
      % According to Marcet and Lorenzoni [1999], if investment >0, then
      % constraints not binding, the same as in general PEA (set mu=0, as constraint are not
      % binding)
      if iv>0;
         % New Capital Stock in t+1
         k(i+1)= (1-delta)*k(i)+iv;
         % Constraint not binding
         mu(i)	= 0;
      else
         % When constraint is binding then calculate capital without
         % investment(negative), consumption as common and mu
         k(i+1)= (1-delta)*k(i);
         % consumption
         c		= exp(a(i))*k(i)^alpha;
         % Lagrange Multiplier
         mu(i)	= c^(-sigma)-lb(i);
      end
   end
   % Our conditional expectation function
   y		= beta*(lb(T1).*(alpha*exp(a(T1)).*k(T1).^(alpha-1)+1-delta)-mu(T1)*(1-delta));
   % As our estimates are well conditioned we stick to OLS, in peaoginit
   % our LS SVD estimate resulted the same b as OLS, so in order to reduce
   % Computational cost we stick to the more simple case OLS
   % Linear regression log(expectation)=const. + bt * log(state variable ...)
   bt		= X(T,:)\log(y);
   % Update paramter b
   b		= gam*bt+(1-gam)*b0;
   % Update coefficient to run loop
   crit	= max(abs(b-b0));
   b0		= b;
   disp(sprintf('Iteration: %d\tConv. crit.: %g',iter,crit))
   iter=iter+1;
end;
% Plot investment
it	= k(T1)-(1-delta)*k(T);
subplot(221)
h=plot(k(T),it,'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

% Plot Histogram of Investment
subplot(222)
hist(it,100);
title('Distribution of investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

% Plot Capital Stock
subplot(223)
h=plot(k(T),k(T1),'.');
set(h,'markersize',0.1,'color','k');
hold on
h1=plot([0:8],(1-delta)*[0:8],'-');
axis([0 8 0 8])
set(h1,'linewidth',1.2,'color','g');
xlabel('k_t','fontname','times','fontsize',12)
title('Capital stock','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

% Plot Lagrange Multiplier
subplot(224)
h=plot(k(T),mu(T),'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('Lagrange multiplier','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

pause
print -depsc2 peairdr.eps
close

% Plot path of investment
T0=500:1000;
subplot(221)
plot(it(T0),'k');
xlabel('Time','fontname','times','fontsize',12)
title('investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

% Plot histogram of Lagrange Multiplier
subplot(222)
plot(mu(T0),'k');
xlabel('Time','fontname','times','fontsize',12)
title('Lagrange multiplier','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

pause
print -depsc2 peairtr.eps
close