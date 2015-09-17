% Diffusion and propulsion - sheared swimmer - Dynamics & EN
% Integrate using Euler-Maruyama with drift correction
% White noise implementation in theta-phi 
% Simplest Explicit Euler used

  clear all;
  RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

  dtmax = 10^-3;                        % max dt
  %nsave = 10^5;                         % save frequency
  sqdt = sqrt(dtmax);            
  nensem = 10^2;                        % number of ensembles

  % Initial orientation and position of COM
  
  theta_0 = pi/4.0                     % Initial theta
  phi_0 = pi/4.0                      % Initial phi

  x_0 = 0.0;  
  y_0 = 0.0;
  z_0 = 0.0;

  % Aspect ratio of particle and run parameters
  
  re  = 10.0;                            % aspect ratio
 aspect_2 = (re^2-1)/(re^2+1);      
  aspect_1 = 1.d0/(re^2+1);

  v0 = 2.0;                              % self propulsion velocity
  B_R = 1.0;                             % strength of Brownian motion_Rotation
  B_T = 1.0;                             % Brownian motion_Translation
  freq = 1.1;                            % Time scaled with inverse frequency 
  Dpar=1.0;                              % D parallel
Tmin=1.0;

  % nstart = floor(Tmin/dt);           % Data collection start time
    
    nstart = 1;
    dt=dtmax;
    Tf = 10.0*Tmin                      
    nT = floor(Tf/dt)  
    C_t0 = zeros(1,nT);
    Cxy = zeros(1,nT);
    cc = zeros(9,nT);
    rr=cc; 
    avgn=zeros(3,nT);
    currn=avgn  ; 
 % Initialize ensemble averaged data
    
   double       nx_0 = sin(theta_0)*cos(phi_0);        
   double       ny_0 = sin(theta_0)*sin(phi_0);         
   double       nz_0 = cos(theta_0);
                    
for B_R=1:10:100;
the=zeros(1,nT);    
for kmk = 1: nensem;   
    
       % Initialize the old values of ensemble for time t = 0
    
       t_o = 0.0;
      % x_o = x_0;
      % y_o = y_0;
      % z_o = z_0;
       phi_o = phi_0;
       the_o = theta_0;
       ncount = 0;
          nx_0 = sin(the_o)*cos(phi_0);         % nx = sin(the)*cos(phi)
          ny_0 = sin(the_o)*sin(phi_0);         % ny = sin(the)*sin(phi)
          nz_0 = cos(the_o);                    % nz = cos(the)

  for  jj = 1: nT;                        % start time stepper for this run 
    
    shear = cos(freq*t_o);            % freq = 1 always in our scalings
    
    nx_o = sin(the_o)*cos(phi_o);         % nx = sin(the)*cos(phi)
    ny_o = sin(the_o)*sin(phi_o);         % ny = sin(the)*sin(phi)
    nz_o = cos(the_o);                    % nz = cos(the)
    
                    % nz = cos(the)
    

     
     % Generate 5 random mumbers
     
         
           
     % Gaussian
     
          
             Ran_the = randn;
             Ran_phi = randn;
     
   % Update positions (x,y,z) - Euler method  
     
   %x_n = x_o + (dt*v0*nx_o) + (dt*shear*y_o) + (sqrt(2.0)*B_T*sqdt*Ran_x);
   %y_n = y_o + (dt*v0*ny_o) + (sqrt(2.0)*B_T*sqdt*Ran_y);
   %z_n = z_o + (dt*v0*nz_o) + (sqrt(2.0)*B_T*sqdt*Ran_z);

   % Update Orientation (the, phi) - Euler method

  the_shear = (dt*0.25)*(aspect_2)*shear*(sin(2.0*the_o)*sin(2.0*phi_o));  
   phi_shear = - (dt*aspect_1)*(shear)*((re*re*sin(phi_o)*sin(phi_o))+(cos(phi_o)*cos(phi_o)));
   
   the_noise = (dt*B_R*B_R*cot(the_o)) + (sqrt(2.0)*B_R*sqdt*Ran_the);   
   phi_noise = +(sqrt(2.0)*B_R*sqdt*Ran_phi)/sin(the_o);  
   
   phi_n = phi_o  + phi_noise; + phi_shear;
   the_n = the_o + the_noise; + the_shear;
   
   t_n = t_o + dt;
 
     %if ((jj > nstart)) %&& (mod(jj,nsave) == 0))         % save
 
                           ncount = ncount + 1;

                           t(ncount) = t_n;       
                           phi(ncount) = phi_n;
                           the(ncount) = the(ncount)+the_n;     
                           
                           %x(ncount) = x_n;
                           %y(ncount) = y_n;
                           %z(ncount) = z_n;
                           
   %  end  
    
     % update old values for next time step
     
     t_o = t_n;
     phi_o = phi_n;
     the_o = the_n;

     %x_o = x_n;
     %y_o = y_n;
     %z_o = z_n;

     

  end      % end of this ensemble run
   
end
the=the/nensem;
%      nx = sin(the).*cos(phi);          % nx = sin(the)*cos(phi)
 %     ny = sin(the).*sin(phi);          % ny = sin(the)*sin(phi)
  %    nz = cos(the);                    % nz = cos(the)
%currn = [sin(the).*cos(phi) sin(the).*sin(phi)  cos(the)];
%currn = [nx; ny; nz];
     
 %c = [nx*nx_0 ; nx*ny_0 ; nx*nz_0; ny*nx_0 ; ny*ny_0 ; ny*nz_0; nz*nx_0 ; nz*ny_0 ; nz*nz_0];
  % cc=cc+c;
  

   % end                                         

      
 %     cc=cc/nensem;
             afname = sprintf('freq%dD_R%densemble.mat',freq,B_R);
              save(afname, 't','the');
end
