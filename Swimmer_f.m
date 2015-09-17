% Stochastic SDE solver for a swimmer in linear oscillatory shear 
% Integrate in the Ito scheme using the  Euler-Maruyama approach
% For shear flow - it is essential to keep the finite aspect ratio if there
% is no brownian motion but we will still keep it as its asymptotically
% correct

% Checked November 12, 2011

clear all;
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

% Set constants 
% Need ideally gamdot to be reasonably small
% Can control either frequency or amplitude

% control time step and ode order

dt=2^-9;
nT=2^18;
T_f=(nT-1)*dt;
sqdt=sqrt(dt);

% number of ensembles

nensem=1000; 

% activity parameters

v0=20.0; 
brow_R = 1.0;
brow_T = 1.0;

% geometry parameters

r_e =20.0 % aspect ratio
aspect=(r_e^2-1.0)/(r_e^2+1.0); %aspect factor

% flow parameters

omega = 1.0; %oscillatory strain rate frequency
epsilon = 0.5; %oscillatory strain rate magnitude
gamdot_0 = 0.0; %constant strain rate

% effective parameters

gamdot = (epsilon*omega+gamdot_0);
gamdot_r1=gamdot/(r_e^2+1.0);
gamdot_r2=gamdot*aspect;

% set initial position and orientation for this bunch of ensembles

x0=0.0;
y0=0.0;
z0=0.0;

theta0=pi/4.0;
phi0=pi/4.0;


% preallocate for ensemble average

MSQX = zeros(1,nT);
MSQY = zeros(1,nT);
MSQZ = zeros(1,nT);

    
 for  kmk = 1:nensem
     
% pre-allocate space for position and orientation   - can remove
% this if needed when memory runs low but at expense of speed

t=zeros(1,nT);

x=zeros(1,nT);
y=zeros(1,nT);
z=zeros(1,nT);

the=zeros(1,nT);
phi=zeros(1,nT);

%  Start run at t=0

t(1)=0.0;
the(1)=theta0;
phi(1)=phi0;

x(1)=x0;
y(1)=y0;
z(1)=z0;

for  jj = 2: nT;
    
    xold=x(jj-1);
    yold=y(jj-1);
    zold=z(jj-1);

    theold=the(jj-1);
    phiold=phi(jj-1);
    
    ct=cos(theold);
    st=sin(theold);
    cp=cos(phiold);
    sp=sin(phiold);
    
    t(jj)=t(jj-1)+dt;
    
    tt=t(jj);
    
    % step forward using Euler-Maruyama
    
    % n_x = sin(the)*cos(phi)
    % n_y = sin(the)*sin(phi)
    % n_z = cos(the)
    % the = azimuthal angle from vorticity axis
    % phi angle made with x-axis of projection in x-y plane

    % dr_i = dt*[vo*(n_i) + gamdot*(ye_y.e_i)] + dW sqrt(dt)
    
    xnew = xold + (dt*v0*(st*cp)) + (brow_T*sqdt*randn) + (dt*gamdot*cos(omega*tt)*yold);
    ynew = yold + (dt*v0*(st*sp)) + (brow_T*sqdt*randn);
    znew = zold + (dt*v0*(ct)) + (brow_T*sqdt*randn);
    
    % Orientation

    % d(the) = dt*[(1/4)*(aspect)*sin(2*the)*sin(2*phi)] + dW 
    % d(phi) = - dt *[gamdot/(r_e^2+1)] * [r_e^2*sin^2(phi) + cos^2(phi)] + dW
    % gamdot = gamdot_0 + e*w*cos(wt)

    phinew = phiold - (dt*(gamdot_r1)*cos(omega*tt)*((r_e^2*sp*sp)+cp*cp)) + (brow_R*sqdt*randn); 
    thenew = theold + (0.25*dt*(gamdot_r2)*cos(omega*tt)*sin(2*theold)*sin(2.0*phiold)) + (brow_R*sqdt*randn); 
    
    % Add new values to the next array position
    
    x(jj)=xnew;
    y(jj)=ynew;
    z(jj)=znew;

    the(jj)=thenew;
    phi(jj)=phinew;
   
end

% end of ensemble_kmk

% Generate trajectories and MSD


  nx=sin(the).*cos(phi);
  ny=sin(the).*sin(phi);
  nz=cos(the);

xx=(x-x0).*(x-x0);
yy=(y-y0).*(y-y0);
zz=(z-z0).*(z-z0);

% Note that the  value of xx(t), yy(t), zz(t) is the MSD evaluated at time t
% This allows us to do both short time and long time MSD's
% over many ensemble averages

    
MSQX = MSQX + xx;
MSQY = MSQY + yy;
MSQZ = MSQZ + zz;

  figure(1)
  plot(t(1:nT),x(1:nT),t(1:nT),y(1:nT), t(1:nT), z(1:nT))
  hold on;

  figure(2)
  plot3(nx(1:nT),ny(1:nT),nz(1:nT))
  hold on;

 end % end of all ensemble runs
 
% average over ensembles

MSQX = MSQX/nensem;
MSQY = MSQY/nensem;
MSQZ = MSQZ/nensem;

 figure(3)
 plot(t(1:nT),MSQX(1:nT))
 title('MSQX')

 
 figure(4);
 plot(t(1:nT),MSQY(1:nT),t(1:nT), MSQZ(1:nT),'--')
 title('MSQY and MSQZ (--)')
 
 
 
 
 

  
