clear

crit	= 1;
tol 	= 1e-6;
gam	= 1;%0.25;
long	= 20000;
init	= 500;
slong	= long+init;
T		= init+1:slong-1;
T1		= init+2:slong;

rw		= 0.7;
sw		= 0.1;
wb		= 0;
beta	= 0.95;
R		= 1/(beta+0.01);
sigma	= 1.5;
ab		= 0;
ass	= 0;

randn('state',1);
e		= sw*randn(slong,1);
w		= zeros(slong,1);
w(1)	= wb+e(1);
for i	= 2:slong;
   w(i)= rw*w(i-1)+(1-rw)*wb+e(i);
end
w=exp(w);

a		= zeros(slong,1);
c		= zeros(slong,1);
lb		= zeros(slong,1);
X		= zeros(slong,6);
a(1)	= 0;
csw	= 0.8;
rt		= 0.2;
sc		= 0.1;
randn('state',1234567890);
ec		= sc*randn(slong,1);
for i=1:slong;
   la		= (a(i));
   X(i,:)= [1 la w(i) la*la w(i)*w(i) la*w(i)];
   c(i)	= rt*a(i)+w(i)+ec(i);
   a1		= R*a(i)+w(i)-c(i);
   if a1>ab;
      a(i+1)=a1;
   else
      a(i+1)= ab;
      c(i)	= R*a(i)+w(i)-ab;
   end
end
lb	= c.^(-sigma);
plot(a(1:200));pause;close
y	= log(beta*R*lb(T1));
b0	= X(T,:)\y

iter=1;
while crit>tol;
   a		= zeros(slong,1);
   c		= zeros(slong,1);
   lb		= zeros(slong,1);
   X		= zeros(slong,length(b0));
   
   a(1)= 0;
   for i=1:slong;
      la		= (a(i));
      X(i,:)= [1 la w(i) la*la w(i)*w(i) la*w(i)];
      lb(i)	= exp(X(i,:)*b0);
      a1		= R*a(i)+w(i)-lb(i)^(-1/sigma);
      if a1>ab;
         a(i+1)=a1;
         c(i)	= lb(i).^(-1./sigma);
      else
         a(i+1)= ab;
         c(i)	= R*a(i)+w(i)-ab;
         lb(i)	= c(i)^(-sigma);
      end
   end
   y	= log(beta*R*lb(T1));
   la= a;
%   ass=mean(a);
   
   b	= X(T,:)\y;
   b	= gam*b+(1-gam)*b0;
   crit=max(abs(b-b0));
   b0=b;
   disp(sprintf('Iteration: %d\tConv. crit.: %g',iter,crit))
   iter=iter+1;
end;

h=plot(R*a(T)+w(T)-ab,c(T),'.');
set(h,'markersize',0.1,'color','k');
xlabel('Cash-on-hand (R a_t+\omega_t-a)','fontname','times','fontsize',12)
title('Consumption','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
pause;
print -depsc2 peabcdr.eps
close

subplot(221)
h=plot(a(T),a(T1),'.');
set(h,'markersize',0.1,'color','k');
xlabel('a_t','fontname','times','fontsize',12)
title('Wealth','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)

subplot(222)
hist(a,100);
title('Distribution of wealth','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
pause;
print -depsc2 peabcdr1.eps
close

