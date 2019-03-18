
function A =Part3() 

    clear all 
    nx = 50;
    ny = 50;
    Ccount = 0; % Collision counter
    count = 0;
    sigma = 10^-2; % box conductivity
    cMap = ones(nx,ny);% The whole region assigned conductivity of 1
    m0 = 9.11e-31; % mass of an electron
    mn = 0.26*m0; % effective electron mass
    kB = 1.38e-23; % boltzman constant
    T = 300; % temperature
    q = 1.602e-19; % electron charge
    Vth = sqrt(2*kB*T/mn);% Thermal Voltage
    xlimit = nx*1e-9; % Region X Limit
    ylimit = ny*1e-9; % Region Y Limit
    n_elec = 1000;
    eleColor = hsv(n_elec);% color palette for electron plotting

    elec = zeros(n_elec, 6);% Electron array for positions and velocities.
    elec_pr = zeros(n_elec, 2); % Previous electron positions
    dt = xlimit/Vth/100; % time step
    totalTime = dt*150; % total length of simulation
    Pscat = 1-exp(-dt/0.2e-12)*2; % scattering probability
    temp= zeros(totalTime/dt,1);% Temperature Array
    timeTrac = zeros(totalTime/dt,1);% Time Array

    elec_density = zeros(100,100); % grid for density map

    % Assign box conductivity valaues to the regions specified by LB and WB.
    LB = 0.4;
    WB = 0.6;
    for m = 1:nx
        for h = 1:ny
            if (m>=LB*nx && m<=WB*nx && h>=WB*ny || m>=LB*nx && m<=WB*nx && h<=LB*ny)
                cMap(m,h) = sigma;
            end
        end
    end

    %% 
    % I will now generate the voltage map of the region using the Finite Difference method
    G = sparse(nx*ny,nx*ny); 
    B= zeros(nx*ny,1);
    for m = 1:nx
        for h = 1:ny
            n = h + (m-1)*ny;

            if m == 1
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = 0.8;

            elseif m == nx
                G(n,:) = 0;
                G(n,n) = 1;

            elseif h == 1
                nxm = h+(m-2)*ny;
                nxp = h+(m)*ny;
                nyp = h+1+(m-1)*ny;

                rxm = (cMap(m,h) + cMap(m-1,h))/2;
                rxp = (cMap(m,h) + cMap(m+1,h))/2;
                ryp = (cMap(m,h) + cMap(m,h+1))/2;

                G(n,n) = -(rxm+rxp+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;

            elseif h == ny
                nxm = h+(m-2)*ny;
                nxp = h+(m)*ny;
                nym = h-1+(m-1)*ny;

                rxm = (cMap(m,h) + cMap(m-1,h))/2;
                rxp = (cMap(m,h) + cMap(m+1,h))/2;
                rym = (cMap(m,h) + cMap(m,h-1))/2;

                G(n,n) = -(rxm+rxp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;

            else 
                nxm = h+(m-2)*ny;
                nxp = h+(m)*ny;
                nym = h-1+(m-1)*ny;
                nyp = h+1+(m-1)*ny;

                rxm = (cMap(m,h) + cMap(m-1,h))/2;
                rxp = (cMap(m,h) + cMap(m+1,h))/2;
                rym = (cMap(m,h) + cMap(m,h-1))/2;
                ryp = (cMap(m,h) + cMap(m,h+1))/2;

                G(n,n) = -(rxm+rxp+ryp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;

            end

        end

    end

    %%
    % The solution is calculated using matrix left division. Then it is mapped from a 1-D column into a 2-D image. 
    V = G\B;

    % Create matrix to map voltage and loop through matrix to assign values
    % from calculated voltage matrix
    VMap = zeros(nx,ny);

    for m = 1:nx
        for h = 1:ny
            n = h + (m-1)*ny;

            VMap(m,h) = V(n);
        end
    end
    VMap = VMap';


    [Ex,Ey] = gradient(-VMap.*1e-9); % Electric Field Gradient Calculation 
    Fx = -Ex*q/5e-8; 
    Fy = -Ey*q/5e-8;
    ax = Fx/mn;
    ay = Fy/mn;
    Vaddx = ax*dt;
    Vaddy =ay*dt;

    figure(9)
    surf(ax)
    title('Acceleration Map')
    xlabel('X nm')
    ylabel('Y nm')
    view(2)


    % Bottle neck co-ordinates
    ybottom = 0.2e-7;
    ytop = 0.3e-7;
    xleft = 0.2e-7;
    xright = 0.3e-7;

    % Initialize electron positions 
    for m=1:n_elec
        for h=1:6
            if(h==1)
                elec(m,h) = xlimit*rand(); % initial x-cordinate 
            elseif(h==2)
                elec(m,h) = ylimit*rand(); % initial y - cordinate
            while((elec(m,1) >= xleft && elec(m,1) <= xright && elec(m,2) <= ybottom) || (elec(m,1) >= xleft && elec(m,1) <= xright && elec(m,2) >= ytop))
                % while electron position is in the square keep generating a new position until its not in the square.
                elec(m,h) = xlimit*rand();
                elec(m,h) = ylimit*rand();
            end
            elseif(h==3)
                elec(m,h) = 2*pi*rand(); % random initial angle
            elseif(h==4)
                elec(m,h) = randn()*Vth; % random initial velocity
            elseif(h==5)
                elec(m,h) = cos(elec(m,3))*elec(m,4)*0.9; % y - velocity component
            else
                elec(m,h) = sin(elec(m,3))*elec(m,4)*0.9; % x - velocity component
            end
        end
    end

    % Elelctron movement
    for g=0:dt:totalTime

        if(Pscat > rand()) % Scattering Check
                vx_new = randn()*Vth/(Ccount*0.5); % velocities decrease over time because of scattering occurences
                vy_new = randn()*Vth/(Ccount*0.5); % velocities decrease over time because of scattering occurences
                elec(m,3) = 2*pi*rand(); % new angle
                v_new = sqrt(vx_new^2+vy_new^2); % Velocity magnitude
                elec(m,4) = v_new;
                elec(m,5) = cos(elec(m,3))*v_new;
                elec(m,5) = sin(elec(m,3))*v_new;
                Ccount =+ 1;
            end

        for m=1:n_elec 
            % if electron is in between the bottle neck boundary condition
            if((elec(m,1) >= xleft && elec(m,1) <= xright) && (elec(m,2) >= ybottom && elec(m,2) <= ytop))
                if(((elec(m,2)+elec(m,6)*dt) >= ytop) || ((elec(m,2)+elec(m,6)*dt) <= ybottom))
                      elec(m,6) = -1*elec(m,6);
                end
            end

           % electron entering box on the left boundary condition
            if(elec(m,1) <= xleft && (elec(m,2) <= ybottom || elec(m,2) >= ytop))
                if((elec(m,1)+elec(m,5)*dt) >= xleft)
                      elec(m,5) = -1*elec(m,5);
                end
            end

            % electron entering box on the right boundary condition
            if(elec(m,1) >= xright && (elec(m,2) <= ybottom || elec(m,2) >= ytop))
                 if((elec(m,1)+elec(m,5)*dt) <= xright)
                      elec(m,5) = -1*elec(m,5);
                end
            end

            if (elec(m,1) >= xlimit) % continous side boundaries
                elec_pr(m,1) = 0;
                elec(m,1) = 0;
            elseif (elec(m,1) <= 0)% continous side boundaries
                elec_pr(m,1) = xlimit;
                elec(m,1) = xlimit;
            end

            if ((elec(m,2) >= ylimit) || (elec(m,2) <= 0)) % reflective top and bottom boundaries
                  elec(m,6) = -1*elec(m,6);
            end 
        end

       % Update previous electron positions
       elec_pr(:,1) = elec(:,1);
       elec_pr(:,2) = elec(:,2);
       % update electron position according to new velocity
       elec(:,1) = elec(:,1) + elec(:,5).*dt;
       elec(:,2) = elec(:,2) + elec(:,6).*dt;

       count = count +1;
       timeTrac(count,1) = g + dt;

       for m = 1:n_elec   
        % find position of the electron in the E-Field
        xpos=round(abs(elec(m,1)*1e9));
        ypos=round(abs(elec(m,2)*1e9));
        b=xpos<50;
        a=xpos~=0;
        c=ypos<50;
        d=ypos~=0;

            if(a && b && c &&d) %% Add the interaction of the E-field 
                elec(m,5) = elec(m,5) + Vaddx(round(abs(elec(m,1)*1e9)),round(abs(elec(m,2)*1e9)))*dt;
                elec(m,6) = elec(m,6) + Vaddy(round(abs(elec(m,1)*1e9)),round(abs(elec(m,2)*1e9)))*dt; 
            end
       end
    end

    % Create electron density map
    for XPos=1:100
        for YPos=1:100
            for q = 1:n_elec
                if((elec(q,1) <= (xlimit*(XPos/100))) && (elec(q,1) > (xlimit*((XPos-1)/100))) && (elec(q,2) <= (ylimit*(YPos/100))) && (elec(q,2) > (ylimit*((YPos-1)/100))))
                    elec_density(XPos,YPos) =+ 1; % tally up eelectrons in the specific position
                end
            end
        end
    end

    figure(10)
    surf(elec_density')
    title('Electron Density Plot')
    xlabel('X 1 square = 0.5nm')
    ylabel('Y 1 square = 0.5nm')
    view(2)

    %% 
    % This code produces the electron trajectory plot for 10 electrons.

    n_elec=10;
    eleColor = hsv(n_elec);
    totalTime=25*dt;
    xpos=0;
    ypos=0;
    a=0;
    b=0;
    c=0;
    d=0;

    for g=0:dt:totalTime
        avg_temp = 0;
        avg_velocity = 0;
        for m=1:n_elec

           if(Pscat > rand()) % Scattering Check
                vx_new = randn()*Vth*2/(Ccount*0.5); % velocities decrease over time because of scattering occurences
                vy_new = randn()*Vth*2/(Ccount*0.5); % velocities decrease over time because of scattering occurences
                elec(m,3) = 2*pi*rand(); % new angle
                v_new = sqrt(vx_new^2+vy_new^2); % Velocity magnitude
                elec(m,4) = v_new;
                elec(m,5) = cos(elec(m,3))*v_new;
                elec(m,5) = sin(elec(m,3))*v_new;
                Ccount =+ 1;
            end

        for m=1:n_elec 
            % if electron is in between the bottle neck boundary condition
            if((elec(m,1) >= xleft && elec(m,1) <= xright) && (elec(m,2) >= ybottom && elec(m,2) <= ytop))
                if(((elec(m,2)+elec(m,6)*dt) >= ytop) || ((elec(m,2)+elec(m,6)*dt) <= ybottom))
                      elec(m,6) = -elec(m,6);
                end
            end

           % electron entering box on the left boundary condition
            if(elec(m,1) <= xleft && (elec(m,2) <= ybottom || elec(m,2) >= ytop))
                if((elec(m,1)+elec(m,5)*dt) >= xleft)
                      elec(m,5) = -elec(m,5);
                end
            end

            % electron entering box on the right boundary condition
            if(elec(m,1) >= xright && (elec(m,2) <= ybottom || elec(m,2) >= ytop))
                 if((elec(m,1)+elec(m,5)*dt) <= xright)
                      elec(m,5) = -elec(m,5);
                end
            end

            if (elec(m,1) >= xlimit) % continous side boundaries
                elec_pr(m,1) = 0;
                elec(m,1) = 0;
            elseif (elec(m,1) <= 0)% continous side boundaries
                elec_pr(m,1) = xlimit;
                elec(m,1) = xlimit;
            end

            if ((elec(m,2) >= ylimit) || (elec(m,2) <= 0)) % reflective top and bottom boundaries
                  elec(m,6) = -elec(m,6);
            end

            %plot the movement of each electron
            if(g~=0)
                 figure(11)
                 plot([elec_pr(m,1),elec(m,1)],[elec_pr(m,2),elec(m,2)],'color',eleColor(m,:))
                 axis([0 xlimit 0 ylimit]);
                 rectangle('Position',[xleft 0 (xright-xleft) ybottom])
                 rectangle('Position',[xleft ytop (xright-xleft) ylimit])
            end
        end

         title('electron trajectory ')
         xlabel('x(nm)')
         ylabel('y(nm)')
         hold on

       % Update previous electron positions
       elec_pr(:,1) = elec(:,1);
       elec_pr(:,2) = elec(:,2);
       % update electron position according to new velocity
       elec(:,1) = elec(:,1) + elec(:,5).*dt;
       elec(:,2) = elec(:,2) + elec(:,6).*dt;

       for m = 1:n_elec   
        % find position of the electron in the E-Field
        xpos=round(abs(elec(m,1)*1e9));
        ypos=round(abs(elec(m,2)*1e9));
        b=xpos<50;
        a=xpos~=0;
        c=ypos<50;
        d=ypos~=0;

            if(a && b && c &&d) %% Add the interaction of the E-field 
                elec(m,5) = elec(m,5) + Vaddx(round(abs(elec(m,1)*1e9)),round(abs(elec(m,2)*1e9)))*0.5e10;
                elec(m,6) = elec(m,6) + Vaddy(round(abs(elec(m,1)*1e9)),round(abs(elec(m,2)*1e9)))*0.5e10; 
            end
       end 
    end

    end
end
