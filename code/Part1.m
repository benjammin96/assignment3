
function D = Part1()
    clear all
    
    m0 = 9.11e-31;
    mn = 0.26*m0;
    kB = 1.38e-23;
    T = 300;
    q = 1.602e-19;
    n = 10^15; 
    xlimit = 200e-9; % region limit
    ylimit = 100e-9; % region limit
    V = 0.1; % voltage applied across semiconductor
    n_elec = 10000; % number of electrons
    vth = sqrt(2*kB*T/mn); % thermal velocity
    eleColor = hsv(n_elec); % color palette
    count = 0;
    dt = xlimit/vth/100; % time step 
    totalTime = dt*500; % total simulation time     
    Pscat = 1-exp(-dt/0.2e-12); % scattering condition
    Ccount = 0;

    elec_pr = zeros(n_elec, 2);% previous electron positions
    elec = zeros(n_elec, 6); % electron information (position, velocity etc)
    
    E = V/xlimit; % Calculated electric field
    F = E*q; % force on electrons
    a = F/mn; % acceleration 
    fprintf("The electric field strength is %d \n",E);
    fprintf("The force on any electron is %d \n",F);
    fprintf("Acceleration of the electron due to the Force %d \n",a);

    temp= zeros(totalTime/dt,1); % temperature array
    timeHolder = zeros(totalTime/dt,1); % time array

% Initialization of electrons
    for i=1:n_elec
        for j=1:6
            if(j==1)
                elec(i,j) = xlimit*rand(); % initial x - cordinate
            elseif(j==2)
                elec(i,j) = ylimit*rand(); % initial y - cordinate
            elseif(j==3)
                elec(i,j) = 2*pi*rand(); % random initial angle
            elseif(j==4)
                elec(i,j) = randn()*vth; % random initial velocity
            elseif(j==5)
                elec(i,j) = elec(i,4)*cos(elec(i,3)); % y - velocity component
            else
                elec(i,j) = elec(i,4)*sin(elec(i,3)); % x - velocity component
            end
        end
    end


    avg_drift = zeros(totalTime/dt,1);
    drift_current = zeros(totalTime/dt,1);
    drift_vx = 0;

% Electron Movement
    for k=0:dt:totalTime
        avg_temp = 0;
        avg_velocity = 0;
        for m=1:n_elec
            
            if(Pscat > rand())  % Scattering Check
                elec(m,3) = 2*pi*rand();
                vx_new = randn()*vth; % velocities decrease over time because of scattering occurences
                vy_new = randn()*vth; % velocities decrease over time because of scattering occurences
                v_new = sqrt(vx_new^2+vy_new^2);
                elec(m,4) = v_new;
                elec(m,5) = cos(elec(m,3))*v_new;
                elec(m,5) = sin(elec(m,3))*v_new;
                 Ccount =+ 1;
            end
       
            if (elec(m,1) >= xlimit) % continous side boundaries
                elec(m,1) = 0;
                elec_pr(m,1) = 0;
                
            elseif (elec(m,1) <= 0) % continous side boundaries
                elec(m,1) = xlimit;
                elec_pr(m,1) = xlimit;
            end
      
            if ((elec(m,2) >= ylimit) || (elec(m,2) <= 0)) % reflective top and bottom boundaries
                elec(m,6) = -elec(m,6);
            end
            
            avg_temp = avg_temp + (elec(m,4)^2)*mn/(2*kB);
            drift_vx = drift_vx + elec(m,5);
        end

       % Update previous electron positions
       elec_pr(:,1) = elec(:,1);
       elec_pr(:,2) = elec(:,2);
       % update electron position according to new velocity
       elec(:,1) = elec(:,1) + elec(:,5).*dt;
       elec(:,2) = elec(:,2) + elec(:,6).*dt;
       
       elec(:,5) = elec(:,5) + a*dt; % add the E-field acceleration

       count = count +1;
       temp(count,1) = avg_temp/n_elec;
       timeHolder(count,1) = k + dt;
       avg_drift(count,1) = drift_vx/n_elec; % average drift velocity
       drift_current(count,1) = q*n*avg_drift(count,1);% average drift current 

    end

       figure(1)
       plot(timeHolder,drift_current)
       title('System Current Density ')
       xlabel('time (s)')
       ylabel('Current (A/m^2)')



    elec_density = zeros(100,100);
    temp_density = zeros(100,100);

   % Create electron density map
    for x_pos=1:100
        for y_pos=1:100
            for q = 1:n_elec
                if((elec(q,1) <= (xlimit*(x_pos/100))) && (elec(q,1) > (xlimit*((x_pos-1)/100))) && (elec(q,2) <= (ylimit*(y_pos/100))) && (elec(q,2) > (ylimit*((y_pos-1)/100))))
                    elec_density(x_pos,y_pos) =+ 1;% tally up eelectrons in the specific position
                    temp_density(x_pos,y_pos) =+ (elec(q,4)^2)*mn/(2*kB); % Calculate the temperature as a function of electrons velocity at that point.
                end
            end
        end
    end



    figure(2)
    surf(elec_density)
    title('Electron Density Map')
    xlabel('x-axis (nm)')
    ylabel('y-axis  (nm)')
    colorbar
    view(2)

    figure(3)
    surf(temp_density)
    title('Temperature Map')
    xlabel('x-axis position (nm)')
    ylabel('y-axis position (nm)')
    colorbar
    view(2)

    %% 

    n_elec=10;
    eleColor = hsv(n_elec);
    totalTime =50*dt;
    
for k=0:dt:totalTime
        avg_temp = 0;
        avg_velocity = 0;
        for m=1:n_elec
            if(Pscat > rand())  % Scattering Check
                elec(m,3) = 2*pi*rand();
                vx_new = randn()*vth; % velocities decrease over time because of scattering occurences
                vy_new = randn()*vth; % velocities decrease over time because of scattering occurences
                v_new = sqrt(vx_new^2+vy_new^2);
                elec(m,4) = v_new;
                elec(m,5) = cos(elec(m,3))*v_new;
                elec(m,5) = sin(elec(m,3))*v_new;
                 Ccount =+ 1;
            end
       
            if (elec(m,1) >= xlimit) % continous side boundaries
                elec(m,1) = 0;
                elec_pr(m,1) = 0;
                
            elseif (elec(m,1) <= 0) % continous side boundaries
                elec(m,1) = xlimit;
                elec_pr(m,1) = xlimit;
            end
      
            if ((elec(m,2) >= ylimit) || (elec(m,2) <= 0)) % reflective top and bottom boundaries
                elec(m,6) = -elec(m,6);
            end
                 if(k~=0)
                 figure(4)
                 plot([elec_pr(m,1),elec(m,1)],[elec_pr(m,2),elec(m,2)],'color',eleColor(m,:))
                 axis([0 xlimit 0 ylimit]);
                 hold on
                 end
        end   
        title('Electron Trajectories')
        xlabel('X nm')
        ylabel('Y nm')
       
        
       % Update previous electron positions
       elec_pr(:,1) = elec(:,1);
       elec_pr(:,2) = elec(:,2);
       % update electron position according to new velocity
       elec(:,1) = elec(:,1) + elec(:,5).*dt;
       elec(:,2) = elec(:,2) + elec(:,6).*dt;
       elec(:,5) = elec(:,5) + a*dt; % add the E-field acceleration

end