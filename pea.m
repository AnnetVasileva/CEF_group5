clear all
%
% Solving the Growth model with PEA
%
%
% 
%
long	= 20000;
init	= 500;
slong	= init+long;
T		= init+1:slong-1;
T1		= init+2:slong;
tol	= 1e-6;
crit	= 1;
gam	= 1;%0.25;

sigma	= 1;
delta	= 0.1;
beta	= 0.95;
alpha	= 0.3;
ab		= 0;
rho	= 0.9;
se		= 0.01;
param	= [ab alpha beta delta rho se sigma long init];
ksy	=(alpha*beta)/(1-beta*(1-delta));
yss	= ksy^(alpha/(1-alpha));
kss	= yss^(1/alpha);
iss	= delta*kss;
css	= yss-iss;
csy	= css/yss;
lss	= css^(-sigma);

randn('state',1);
e		= se*randn(slong,1);
a		= zeros(slong,1);
a(1)	= ab+e(1);
for i	= 2:slong;
   a(i)=rho*a(i-1)+(1-rho)*ab+e(i);
end

b0		= peaoginit(e,param);

%
% Main Loop
%
iter	= 1;
while crit>tol;
   %
   % Simulated path
   %
   k		= zeros(slong+1,1);
   lb		= zeros(slong,1);
   X		= zeros(slong,length(b0));
   k(1)	= kss;
   for i	= 1:slong;
      X(i,:)= [1 log(k(i)) a(i) log(k(i))*log(k(i)) a(i)*a(i) log(k(i))*a(i)];
      lb(i)	= exp(X(i,:)*b0);
      k(i+1)=exp(a(i))*k(i)^alpha+(1-delta)*k(i)-lb(i)^(-1/sigma);
   end
   y		= beta*lb(T1).*(alpha*exp(a(T1)).*k(T1).^(alpha-1)+1-delta);
   bt		= X(T,:)\log(y);
   b		= gam*bt+(1-gam)*b0;
   crit	= max(abs(b-b0));
   b0		= b;
   disp(b0(:)')
   disp(sprintf('Iteration: %d\tConv. crit.: %g',iter,crit))
   iter=iter+1;
end;

h=plot(k(T),k(T1),'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
ylabel('k_{t+1}','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
pause
%print -depsc2 peaogdrk.eps
close