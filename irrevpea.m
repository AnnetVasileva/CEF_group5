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
gam	= 0.25;

sigma	= 1;
delta	= 0.1;
beta	= 0.95;
alpha	= 0.3;
ab		= 0;
rho	= 0.8;
se		= 0.125;
kss	= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
css	= kss^alpha-delta*kss;
lss	= css^(-sigma);
ysk	= (1-beta*(1-delta))/(alpha*beta);
csy	= 1-delta/ysk;
%
% Initial guess
%
randn('state',1);
e		= se*randn(slong,1);
a		= zeros(slong,1);
a(1)	= ab+e(1);
for i	= 2:slong;
   a(i)=rho*a(i-1)+(1-rho)*ab+e(i);
end
param	= [ab alpha beta delta rho se sigma long init];

b0		= peaoginit(e,param);

%
% Main Loop
%
iter		= 1;
while crit>tol;
   %
   % Simulated path
   %
   k		= zeros(slong+1,1);
   lb		= zeros(slong,1);
   mu		= zeros(slong,1);
   X		= zeros(slong,length(b0));
   k(1)	= kss;
   for i	= 1:slong;
      X(i,:)= [1 log(k(i)) a(i) log(k(i))*log(k(i)) a(i)*a(i) log(k(i))*a(i)];
      lb(i)	= exp(X(i,:)*b0);
      iv		= exp(a(i))*k(i)^alpha-lb(i)^(-1/sigma);
      if iv>0;
         k(i+1)= (1-delta)*k(i)+iv;
         mu(i)	= 0;
      else
         k(i+1)= (1-delta)*k(i);
         c		= exp(a(i))*k(i)^alpha;
         mu(i)	= c^(-sigma)-lb(i);
      end
   end
   y		= beta*(lb(T1).*(alpha*exp(a(T1)).*k(T1).^(alpha-1)+1-delta)-mu(T1)*(1-delta));
   bt		= X(T,:)\log(y);
   b		= gam*bt+(1-gam)*b0;
   crit	= max(abs(b-b0));
   b0		= b;
   disp(sprintf('Iteration: %d\tConv. crit.: %g',iter,crit))
   iter=iter+1;
end;

it	= k(T1)-(1-delta)*k(T);
subplot(221)
h=plot(k(T),it,'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(222)
hist(it,100);
title('Distribution of investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

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

subplot(224)
h=plot(k(T),mu(T),'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('Lagrange multiplier','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

pause
print -depsc2 peairdr.eps
close

T0=500:1000;
subplot(221)
plot(it(T0),'k');
xlabel('Time','fontname','times','fontsize',12)
title('investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(222)
plot(mu(T0),'k');
xlabel('Time','fontname','times','fontsize',12)
title('Lagrange multiplier','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

pause
print -depsc2 peairtr.eps
close