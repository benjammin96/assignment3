%% Part 1 
% To begin the assignment , the first task was a Monte - Carlo simulation
% of electrons moving in an empty region in the influence of a electric
% field. 
%
% The magnitude of electric field present in the region was 50,000 V/m
% given by:
%
% $$E=\frac{V}{\Delta X}$$
%
% Where V -0.1 and Delta X =200nm.
%
% Each electron felt a force of 8.1E-14 N due to this electric field.
%
% $$F=q\times E$$
%
% and the electrons experienced an acceleration of 3.82E17 m/s^2.
%
% $$a= \frac{F}{m}$$
%
% The average current density of the system wasfound using:
%
% $$J = q\times n\times \bar{v}$$
%
% n: electron density 1*10^15 cm^2 , q: elementary charge 1.602E-19C, vbar:
% average drift velocity in the system.
%
% The electron trajectories, the current density plot, electron
% density map and temperature map can be found in this section.
% 
close all 
clear all
%set(groot,'defaultFigureVisible','off')% Figures off
set(groot,'defaultFigureVisible','on') % Figures on
Part1()
%% Part 1 Summary 
% The electric field strength is 50,000 V/m of a 200nm region length. The
% force on any electron is 8.1E-14 N and the acceleration of any electron
% is 3.382E17 m/s^2. Due to this acceleration it is observed that the
% electrons paths are curving. 
%
% The systems electric field is constantly accelerating the electrons and
% thus increasing the current density. The electrons are distributed
% uniformly throughout the region as seen in the electron density and
% temperature plots. 

%% Part 2 
% In this part we computed the voltage potential of the entire region using
% the Finite Difference method. The boundary conditions were set as
% follows.
%
% V =0.8volt at x=0nm and V=0 at x=50nm. The top and bottom boundaries
% were set such that their first derivative was 0.
% 
% The conductivity of the region was initialized as 1 except for the area
% which the boxes contained having a conducvity of 0.01.

Part2();

%% Part 2 Summary 
% The voltage changes rapidly near the lower conductivity boxes to satisfy
% the boundary conditions. In addition, we see that the electric field is
% quite strong between the two boxes and near the corners. 


%% Part 3 
% The goal of part 3 was to incorporate the electric field results of part 2 as an input
% to the monte carlo simulation of electron movements from part 1. 
%
% The electric field influences the electron movement because the forces
% acting on the electrons are described as:
%
% Force:
%
% $F=E\times q$ 
%
% Acceleration:
%
% $a = \frac{F}{m}$
%
% Additive Speed:
%
% $V_{addtivie}=a\times dt$
% 
% 
% A plot of electron density, electron trajectories and the
% acceleration due to the electric field result are displayed.

Part3()

%% Part 3 Summary 
% Observing the electron density plot we can see that there are no
% electrons in the boxes as we expect due to to their opaque boundaries. In
% addition we can see an even distribution of electrons throughout the
% plot.  This makes sense because the increased acceleration due to the
% E-field symmetrical about a y -axis in the middle of the plot. 
%
% To increase the simulations accuracy we could add more electrons to the
% simulation and increase the mesh density during the finite difference
% calculation of the voltage solution.
% 
%