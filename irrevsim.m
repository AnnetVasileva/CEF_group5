clear all
%
% Solving the Growth model with PEA
%
%
% 
%
long	= 200;
slong	= long;
T		= 1:slong-1;
T1		= 2:slong;

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
e		= -1.5*se*ones(slong,1);
a		= zeros(slong,1);
a(1)	= ab+e(1);
for i	= 2:slong;
   a(i)=rho*a(i-1)+(1-rho)*ab+e(i);
end

b0 =[
   0.35579546312586
  -0.32885591385471
  -0.71819673342925
  -0.12009049070903
  -0.21676093877187
   0.31262562596598];
%
% Main Loop
%
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
   
   it	= k(T1)-(1-delta)*k(T);
subplot(221)
h=plot(k(T),it);
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(222)
hist(it,100);
title('Distribution of investment','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(223)
h=plot(k(T),k(T1));
set(h,'markersize',0.1,'color','k');
%hold on
%h1=plot([0:8],(1-delta)*[0:8],'-');
%axis([0 8 0 8])
%set(h1,'linewidth',1.2,'color','g');
xlabel('k_t','fontname','times','fontsize',12)
title('Capital stock','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(224)
h=plot(mu(T));
set(h,'markersize',0.1,'color','k');
xlabel('k_t','fontname','times','fontsize',12)
title('Lagrange multiplier','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

pause
%print -depsc2 peairdr.eps
close

