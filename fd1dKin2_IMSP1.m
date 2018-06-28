function fd1dKin2_IMSP1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Discussion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fd1dKin2_IMSP1.m   1D finite-difference method for IMSP first order scheme  
% applied to the predator-prey system with Kinetics 2.
% (C) 2005 Marcus R. Garvie. See COPYRIGHT.TXT for details.
%
% Values to be assigned for reproducing experiments in 
% F.Diele, M.R. Garvie, C.Trenchea, 
% 'Analysis of first-order in time implicit-symplectic scheme
% for predator-prey systems', 2016 submitted
%
% alpha=1.5; beta=1; gamma=5; Du=1; Dv=1;
% a=0; b=1; 
% u0='0.2'; v0='0.0328';
% Unstable steady state 
% ustar= -1/gamma*log((alpha-1)/alpha);
% vstar= ustar*(1-ustar)/(1-exp(-gamma*ustar));


% User inputs of parameters
alpha = input('Enter parameter alpha   ');
beta = input('Enter parameter beta   ');
gamma = input('Enter parameter gamma   ');
Du = input('Enter diffusion coefficient Du   ');
Dv = input('Enter diffusion coefficient Dv   ');

a = input('Enter a in [a,b]   ');
b = input('Enter b in [a,b]  ');
T = input('Enter maximum time T   ');


% User inputs of spatial and temporal stepsizes
h = input('Enter space-step h   ');
delt = input('Enter time-step Delta t   ');

% User inputs of initial data
u0 = input('Enter initial data function u0(x)   ','s'); % see notes
v0 = input('Enter initial data function v0(x)   ','s'); %   in text


% Values to be assigned for reproducing experiments in 
% F.Diele, M.R. Garvie, C.Trenchea, 
% 'Analysis of first-order in time implicit-symplectic scheme
% for predator-prey systems', 2016 submitted
%
% alpha=1.5; beta=1; gamma=5; Du=1; Dv=1;
% a=0; b=1; 
% u0='0.2'; v0='0.0328';
% Unstable steady state 
% ustar= -1/gamma*log((alpha-1)/alpha);
% vstar= ustar*(1-ustar)/(1-exp(-gamma*ustar));


% Calculate some constants
mu=1/(h^2);
J=round((b-a)/h);
n = J+1;   % no. of nodes (d.f.) for each dependent variable                        
N=round(T/delt);


% Initialization
u=zeros(n,1); v=zeros(n,1); F=zeros(n,1); G=zeros(n,1); indexI=zeros(n,1);
y1=zeros(n,1); y2=zeros(n,1); z1=zeros(n,1); z2=zeros(n,1); hhat=zeros(n,1);
x=zeros(n,1); B1=sparse(n,n); B2=sparse(n,n); L=sparse(n,n); Lower1=sparse(n,n);
Lower2=sparse(n,n); Upper1=sparse(n,n); Upper2=sparse(n,n);
% Assign initial data
indexI=[1:n]'; 
x=(indexI-1)*h+a;   % vector of x values on grid
u = eval(u0).*ones(n,1); v = eval(v0).*ones(n,1);

% Construct matrix L (without 1/h^2 factor)
L=sparse(1,2,2,n,n);
L=L+sparse(n,n-1,2,n,n);
L=L+sparse(2:n-1,3:n,1,n,n);
L=L+sparse(2:n-1,1:n-2,1,n,n);
L=L+sparse(1:n,1:n,-2,n,n);
% Construct matrices B1 & B2
B1=sparse(1:n,1:n,1,n,n) - Du*mu*delt*L;
B2=sparse(1:n,1:n,1,n,n) - Dv*mu*delt*L;
% Perform the LU factorisation of B1 and B2 
[Lower1,Upper1]=lu(B1);
[Lower2,Upper2]=lu(B2);

 t(1)=0; UU=u; VV=v;
% Time-stepping procedure
for nt=1:N
    % Evaluate modified functional response
    hhat = ones(n,1)-exp(-gamma*abs(u));
    % Update right-hand-side of linear system 
    G = beta*alpha*hhat - beta;
    v = v./(1-delt*G);
    F = u - u.*u - v.*hhat;
    y1 = u + delt*F;
    y2 = v;
    % Forward substitution to solve Lower1*z1=y1 for z1
    z1 = Lower1\y1;
    % Back-substitution to solve Upper1*u=z1 for u
    u = Upper1\z1;
    % Forward substitution to solve Lower2*z2=y2 for z2
    z2 = Lower2\y2;
    % Back-substitution to solve Upper2*v=z2 for v
    v = Upper2\z2; 
    
    
    t(nt+1)=nt*delt;
    UU=[UU,u];
    VV=[VV,v];
end

% Plot solution at time level T=N*delt
% figure
% plot(x,u,'k'); hold on; plot(x,v,'k-.')


%Plot prey densities u in [0,T]
figure
surf(x,t,UU');
title('u(x,t)');
xlabel('Distance x');
ylabel('Time t');
hold on

%Plot predator densities v in [0,T]
figure
surf(x,t,VV');
title('v(x,t)');
xlabel('Distance x');
ylabel('Time t');
hold on
 

