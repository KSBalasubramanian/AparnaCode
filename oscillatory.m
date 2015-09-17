% Sheared Swimmer simulation via Euler-Maruyama
 
  clear all;
%  RandStream.setDefaultStream ... 
  %   (RandStream('mt19937ar','seed',sum(100*clock))); % Seed random legacy code for older Matlab systems

%Cluster specific code-------
%ms.UseParallel = 'always'
%matlabpool close force local
%matlabpool open 8 
 
              %c=clock;          
              %afname =sprintf('OS_%d_%d_%d',c(4),c(3),c(2)); % File name 
              %mkdir(afname);
              %cd(afname);
%Cluster specific code-------


  
  % System params
  
  h  = 10.0;                            % aspect ratio
  aspect_2 = (h^2-1)/(h^2+1);      
  aspect_1 = 1.d0/(h^2+1);
  v0 = 10.0;                              % self propulsion velocity
  
  freq = 0.05;                            % shear frequency
  amp=1.0;                               % shear amplitude
  Dpar=1.0;                              % D Translation parallel
  B_R=1.0;                               % D Rotation

 % Simulation params
                           
  dt=10^-2;                             % timestep
  sqdt=sqrt(dt);   
  nensem = 10;                        % number of ensembles
  
    
   nT = 100/dt;                         %number of runs : 50 D, looks at props at infinity  
  finr=zeros(3,nensem); 
% Initial orientation and position of COM  
  
  theta_0 = 3*pi/8.0;                     % Initial theta
  phi_0 = pi/4.0;                      % Initial phi

  x_0 = 0.0;  
  y_0 = 0.0;
  z_0 = 0.0;
%Randomize if need be

% Initialize ensemble averaged data
    
          nx_0 = sin(theta_0)*cos(phi_0);        
          ny_0 = sin(theta_0)*sin(phi_0);         
          nz_0 = cos(theta_0);                    

for kmk = 1: nensem;   % Main loop : rewite to use multi cores.
    
       % Initialize the old values of ensemble for time t = 0
        t_o = 0.0;
        r_n = [0,0,0];
        r_o=[0,0,0];
 
       phi_o = phi_0;              %rewrite to random
       the_o = theta_0;
       
  for  jj = 1: nT;                        % start time stepper for this run 
    
    shear = amp * sin(freq*t_o);            
    
         
   % Update Orientation (the, phi) - Euler method
    Ran_the = randn;                % Gaussian noise for theta
    Ran_phi = randn;
     
    the_shear = (dt*0.25)*(aspect_2)*shear*(sin(2.0*the_o)*sin(2.0*phi_o));  
    phi_shear = - (dt*aspect_1)*(shear)*((h*h*sin(phi_o)*sin(phi_o))+(cos(phi_o)*cos(phi_o)));
   
    the_noise = (dt*B_R*B_R*cot(the_o)) + (sqrt(2.0)*B_R*sqdt*Ran_the);   
    phi_noise = +(sqrt(2.0)*B_R*sqdt*Ran_phi)/sin(the_o);  
   
    phi_n = phi_o  + phi_noise + phi_shear;
    the_n = the_o + the_noise + the_shear;
  
           
   % Update positions (x,y,z) - Euler method
    randpar=randn;
    randperp=randn; % two gaussian noise for r 
    n = [sin(the_o)*cos(phi_o),sin(the_o)*sin(phi_o),cos(the_o)];   %n vector
  
    r_noise=[0,0,0];                      %r noise vector
    perp1 = r_o - dot(r_o,n);              %perp 1 vector
    perp2 = cross(n,perp1);                %perp 2 vector

    n_par =  +(sqrt(2.0))*Dpar*sqdt*randpar; %parallel noise 
    n_perp= +(sqrt(2.0))*0.5*Dpar*sqdt*randperp; %perpendicular noise D_perp is 0.5 D_par

    r_noise= n_par*n + n_perp*(perp1+perp2)  ;
    
    r_n = r_o + dt*v0*n +r_noise+dt*shear*dot(r_o,[0,1,0])*[1,0,0];  %shear term need y_0
   
      
    t_n = t_o + dt;
 
     % update old values for next time step
     
     t_o = t_n;
     phi_o = phi_n;
     the_o = the_n;
     r_o=r_n; 
     t(jj)=t_n;
          
  end      % end of this ensemble run
         r_o
         r_n 
         finr(:,kmk)=r_o
         
end    % finished all ensembles                                     

              c=clock;           
              afname =sprintf('freq%dv%d_%d_%d_%d',freq,v0,c(4),c(3),c(2)); % File name OscSwimmer date time
              save(afname, 'finr'); 
%To save other variables just add them in save save(afname,'finr','var2','var3')
