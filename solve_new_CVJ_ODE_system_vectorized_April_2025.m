
% Program to solve my new generalized ODE. Clifford V. Johnson 13th-16th April 2025
% Vectorized, and with Jacobians. 
% I use chebfuns to smoothly handle some functions but unneccessary in this version.
% Will need the chebfun library, or just take them out.

clear all;
close all;

global xL xR h  r1cheb1 r1cheb2 r2cheb1 r2cheb2 scheb1 scheb2  E1 E2 

 
N=20000;
MaxN=100000;

h = 1;

E1 = 1;

E2 = 1.01;




 
xL = -40; xR = 40;
x = linspace(xL,xR,N); 


%% Set up asymptotic behaviour, to be used as initial guess. 


disp('solving CVJs new ODE system');


disp('constructing guess functions');

%The example here is Airy model, so u(x) = -x exactly.


% first for R1
r1 = -1/2 * 1./ sqrt(-x - E1) ;
r1diff = -(-x - E1) .^ (-0.3e1 / 0.2e1) / 0.4e1;

% singularity in middle for guess seems (later) to be not a problem


r1fun = @(y) interp1(x,r1,y,'spline');
r1difffun = @(y) interp1(x,r1diff,y,'spline');
r1cheb1 = chebfun(r1fun,[xL,xR],'splitting','on');
r1cheb2= chebfun(r1difffun,[xL,xR],'splitting','on');

% I convert them to chebfun for smoothness if I do multiple
% differentiations (on more complicated u(x) examples later)
% I can use these functions as input in the initial guess later


% second for R2 
r2 = -1/2 * 1./ sqrt(-x - E2) ;
r2diff = -(-x - E2) .^ (-0.3e1 / 0.2e1) / 0.4e1  ;

r2fun = @(y) interp1(x,r2,y,'spline');
r2difffun = @(y) interp1(x,r2diff,y,'spline');
r2cheb1 = chebfun(r2fun,[xL,xR],'splitting','on');
r2cheb2= chebfun(r2difffun,[xL,xR],'splitting','on');


% then for S

% (This is a guess that works well enough for large nagative x. Adapted from analytic result.)

s = -1/2 * 1./ (-x - E1).^(1/4)./(-x - E2).^(1/4) ;

sdiff =  ((2 * x + E1 + E2) .* (-x - E1) .^ (-0.5e1 / 0.4e1) .* (-x - E2) .^ (-0.5e1 / 0.4e1)) / 0.8e1;

sfun = @(y) interp1(x,s,y,'spline');
sdifffun = @(y) interp1(x,sdiff,y,'spline');
scheb1 = chebfun(sfun,[xL,xR],'splitting','on');
scheb2 = chebfun(sdifffun,[xL,xR],'splitting','on');

% Plots in case I want to look at it.
% guessfuncS = figure;
% figure(guessfuncS)
% hold on
% plot(x,scheb1(x)),shg
% plot(x,scheb2(x)),shg
% hold off

%% Solve the equation
disp('Now beginning main solving routine');

solinit = bvpinit(x, @guess);

options = bvpset('FJacobian',@jac,'Stats','on','RelTol',1e-4,'abstol',1e-5,'Vectorized','on','Nmax',MaxN);

tic 

sol = bvp4c(@bvpfcn, @bcfcn, solinit, options);

toc

R1 = deval(sol,x,1);  % the R for E1 solution
R2 = deval(sol,x,4);  % the R for E2 solution
S  = deval(sol,x,7);  % the S solution


%% Finally, look at solution

 % Test it against exact result for Airy model.

  psi1 = h^(-2/3) * airy(-h^(-2/3)*(E1+x));
  psi2 = h^(-2/3) * airy(-h^(-2/3)*(E2+x));

  %temp grid for plotting purposes.


%Overall normalization can be off, so normalize to exact result before comparing.

norm = max(max(psi1.*psi2/imag(S)))

% If solution is mostly good but has bad parts to the extreme right, this will fail. 
% Instead can pick a specific point on grid by hand
% norm = max(max(psi1(4000)*psi2(4000)/imag(S(4000))))

disp('Displaying result');

result = figure;
figure(result)

hold on 
  
plot(x,norm*imag(S),'-k','LineWidth',2),shg
plot(x,norm*imag(R1),'-g','LineWidth',2),shg
plot(x,norm*imag(R2),'-b','LineWidth',2),shg

%temporarily make grid y for a less fine grid
y = linspace(xL,xR,3000); 

psi1 = h^(-2/3) * airy(-h^(-2/3)*(E1+y));
psi2 = h^(-2/3) * airy(-h^(-2/3)*(E2+y));

plot(y,psi1.*psi2,'+r','LineWidth',1),shg

%if curious to see what the real parts are doing
%plot(x,norm*real(S),'--k','LineWidth',2),shg
%plot(x,norm*real(R2),'--b','LineWidth',2),shg
%plot(x,norm*real(R1),'--g','LineWidth',2),shg

% plot S again just to get it on top of the crosses, but keep earlier plot for legend order!
plot(x,norm*imag(S),'-k','LineWidth',2),shg 

% mess around by hand to get nice plot range, depending upon energy choices.
axis([-10 10 min(imag(S))-0.02 max(imag(S))+0.02]);
axis([-10 10 min(norm*imag(R2))-0.015 max(norm*imag(R2))+0.01]);
%axis([-10 10 min(imag(R2))-2 max(imag(R2))+2]);
hold off

xlabel('$x$','Fontsize',25,'Interpreter','Latex'); ylabel(' ','Fontsize',25,'rotation',0,'Interpreter','Latex');
legend('${\rm Im}[{\widehat S}(E,E^\prime,x)]$','${\rm Im}[{\widehat R}(E,x)]$','${\rm Im}[{\widehat R}(E^\prime,x)]$','$\psi(E,x)\psi(E^\prime,x)$','Fontsize',20,'Interpreter','Latex')


%% rejoice!
load handel.mat;
sound(y);

%% What follows are the needed input functions.

%---------------------------------

function dydz = bvpfcn(x,y) % The equation to solve
global h   E1 E2  

 dydz = [y(2,:)
         y(3,:)
         (4*(-x - E1).*y(2,:)+2.*(-1).*y(1,:))/h^2
         y(5,:)
         y(6,:)
         (4*(-x - E2).*y(5,:)+2.*(-1).*y(4,:))/h^2
         y(8,:)
         y(9,:)
         (2*(-x - (E1+E2)/2).*y(8,:)+2*(-1).*y(7,:))/h^2+(1./y(7,:)/h^2).*((-x-E1).*y(1,:).*y(5,:)+(-x-E2).*y(4,:).*y(2,:))];

end

%---------------------------------

function dfdy = jac(z,y) % The Jacobians
global h E1 E2
 
dfdy = [0 1 0 0 0 0 0 0 0
        0 0 1 0 0 0 0 0 0
        1/h^2*(-2) 1/h^2*(4*(-z-E1)) 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0
        0 0 0 0 0 1 0 0 0
        0 0 0 1/h^2*(-2) 1/h^2*(4*(-z-E2)) 0 0 0 0
        0 0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 0 1
        1/y(7)/h^2*(-z-E1)*y(5) 1/y(7)/h^2*(-z-E2)*y(4) 0 1/y(7)/h^2*(-z-E2)*y(2) 1/y(7)/h^2*(-z-E1)*y(1) 0 -2/h^2-1/y(7)^2/h^2 *((-z-E1)*y(1)*y(5)+(-z-E2)*y(4)*y(2)) 2*(-z-(E1+E2)/2)/h^2 0
        ];
end

%---------------------------------

function res = bcfcn(ya,yb) % boundary conditions

global xL xR   h E1 E2 % scheb1 scheb2
res = [ya(1)+1/(-xL-E1)^(0.5)/2
       ya(2)-1/(-xL-E1)^(1.5)/4
       yb(1)+1/(-xR-E1)^(0.5)/2
       ya(4)+1/(-xL-E2)^(0.5)/2
       ya(5)-1/(-xL-E2)^(1.5)/4
       yb(4)+1/(-xR-E2)^(0.5)/2
       ya(7)-(-1/2 * 1./ (-xL - E1).^(1/4)./(-xL - E2).^(1/4))
       ya(8)-(((2 * xL + E1 + E2) .* (-xL - E1) .^ (-0.5e1 / 0.4e1) .* (-xL - E2) .^ (-0.5e1 / 0.4e1)) / 0.8e1)
       yb(7)-(-1/2 * 1./ (-xR - E1).^(1/4)./(-xR - E2).^(1/4))
       ];
end

%---------------------------------

function the_guess = guess(z) % initial guess for y and y' etc

global  scheb1 scheb2 r1cheb1 r1cheb2 r2cheb1 r2cheb2

the_guess = [
    r1cheb1(z)
    r1cheb2(z)
     0
    r2cheb1(z)
    r2cheb2(z)
     0
    scheb1(z)
    scheb2(z)
     0];
end

%---------------------------------