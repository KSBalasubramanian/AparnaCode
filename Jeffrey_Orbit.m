% Jeffrey's orbit - vizualization - NO NOISE
% Arvind Gopinath October, 2011

clear all;
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));
 
% Parameters 

dt = 10^-5;                      % time step
nT = 4*10^6;                   % number of steps 
nsave = 10^3;                 % save frequency

 re = 10;                           % aspect ratio

 freq =  0.0;                     % Oscillatory strain rate frequency
 amp = 5.0;                      % Oscillatory strain rate magnitude 

theta_0 = pi/3.0;             % Initial theta value
phi_0 = pi/4.0;                % Initial phi value

% Start by initialing the code

t_o = 0.0;                 
phi_o = phi_0;
the_o = theta_0;

t(1) = 0.0;
phi(1) = phi_0;
the(1) = theta_0;

ncount = 1;

% Time stepper starts here - using a simple first order Euler code. If
% needed we can convert to second order code bu this is easier for now to 
% play with

for  jj = 2: nT;
    
    shear = amp*cos(freq*t_o);                   % Shear
    
    nxo = sin(the_o)*cos(phi_o);                % x component of director
    nyo = sin(the_o)*sin(phi_o);                 % y component of director
    nzo = cos(the_o);                                   % z component of director 
  
    phi_shear = - (dt/(re^2+1))*(shear)*((re*re*sin(phi_o)*sin(phi_o))+cos(phi_o)*cos(phi_o));         % increment in phi
    
    phi_n = phi_o + phi_shear; 
    
    the_shear = (dt*0.25)*((re^2-1)/(re^2+1))*shear*(sin(2.0*the_o)*sin(2.0*phi_o));                         % increment in theta
          
    the_n = the_o + the_shear; 

    t_n = t_o + dt;                                                                                                                                    % increment in time
  
          if(mod(jj,nsave)==0)

          ncount = ncount + 1;
              t(ncount) = t_n;
              phi(ncount) = phi_n;
              the(ncount) = the_n;
 
      end
                
t_o = t_n;
phi_o = phi_n;
the_o = the_n;

end 

nx = sin(the).*cos(phi);
ny = sin(the).*sin(phi);
nz = cos(the);

figure(1);
plot3(nx,ny,nz);

figure(2);
plot(t,the,t,phi,'--');


