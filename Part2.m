
function C=Part2()
    clear all

    nx = 50; % mesh density
    ny = 50; % mesh density
    LB = 0.4; % Boxes dimensions
    WB = 0.6; % Boxes dimensions
    sigma = 0.01;
    B = zeros(nx*ny,1);
    G = sparse(nx*ny,nx*ny);
    V0=0.8; % Voltage of left side boundary

    % Create conductivity map and initialize the lower conductivity box
    % regions.
    cMap = ones(nx,ny);
    for m = 1:nx
        for h = 1:ny
            if (m>=LB*nx && m<=WB*nx && h>=WB*ny || m>=LB*nx && m<=WB*nx && h<=LB*ny)
                cMap(m,h) = sigma;
            end
        end
    end
% finite Difference method 
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1)*ny;

            if i == 1
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = V0;

            elseif i == nx
                G(n,:) = 0;
                G(n,n) = 1;

            elseif j == 1
                nxm = j+(i-2)*ny;
                nxp = j+(i)*ny;
                nyp = j+1+(i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                ryp = (cMap(i,j) + cMap(i,j+1))/2;

                G(n,n) = -(rxm+rxp+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;

            elseif j == ny
                nxm = j+(i-2)*ny;
                nxp = j+(i)*ny;
                nym = j-1+(i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                rym = (cMap(i,j) + cMap(i,j-1))/2;

                G(n,n) = -(rxm+rxp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;

            else 
                nxm = j+(i-2)*ny;
                nxp = j+(i)*ny;
                nym = j-1+(i-1)*ny;
                nyp = j+1+(i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                rym = (cMap(i,j) + cMap(i,j-1))/2;
                ryp = (cMap(i,j) + cMap(i,j+1))/2;

                G(n,n) = -(rxm+rxp+ryp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;
            end
        end
    end

   % Calculate voltage map using matrix left division
    V = G\B;

    % Map voltage solution into a region map
    VMap = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1)*ny;
            VMap(i,j) = V(n);
        end
    end
    VMap = VMap';
    
        figure(5)
        contourf(VMap,50)
        colorbar
        title('Voltage Distribution')
        xlabel('Y')
        ylabel('X')
        
        figure(7)
        surf(VMap)
        colorbar
        title('Voltage Distribution')
        xlabel('Y')
        ylabel('X')
        view([45 45])
    
        % Calculate electric field gradient. 
        [Ex,Ey] = gradient(VMap);
        Ex=Ex;
        Ey=Ey;

        % Electric field gradient map
        figure(8)
        quiver(Ex,Ey)
        title ('Electric Field Vectors')
        xlabel('Y')
        ylabel('X')
        ylim([0 50])
        xlim([0 50])
end
