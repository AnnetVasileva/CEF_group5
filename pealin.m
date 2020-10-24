clear all
%
% Solving the Growth model with PEA
%
%
% 
%
long	= 100000;
init	= 500;
slong	= init+long;
T		= init+1:slong-1;
T1		= init+2:slong;
tol	= 1e-6;
crit	= 1;
gam	= 1;%0.25;

c0		= 0.25;
c1		= 0.40;
ab		= 0;
rho	= 0.9;
se		= 0.01;

randn('state',1);
e		= randn(slong,1);
e		= se*(e-mean(e))/std(e);
a		= zeros(slong,1);
a(1)	= ab+e(1);
for i	= 2:slong;
   a(i)=rho*a(i-1)+(1-rho)*ab+e(i);
end

b0		= 0.1*ones(2,1);
%
% Main Loop
%
iter	= 1;
while crit>tol;
   %
   % Simulated path
   %
   y		= zeros(slong,1);
   X		= zeros(slong,length(b0));
%   for i	= 1:slong;
      X= [ones(slong,1) a];% a(i)*a(i)];
      y= X*b0;
%   end
   R		= c0*y(T1)+c1*a(T);
   bt		= X(T,:)\R;
   b		= gam*bt+(1-gam)*b0;
   crit	= max(abs(b-b0));
   b0		= b;
   disp(b0(:)')
   disp(sprintf('Iteration: %d\tConv. crit.: %g',iter,crit))
   iter=iter+1;
end;
break
h=plot(k(T),k(T1),'.');
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
ylabel('k_{t+1}','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
pause
%print -depsc2 peaogdrk.eps
close