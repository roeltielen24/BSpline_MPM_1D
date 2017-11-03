%% Oedometer_Bspline Version 19-10-2016
%  Roel Tielen (Based on code of Lisa Wobbes), TU Delft
%  This file provides input for MPM to compute numerical solution, 
%  illustrates results and provides the RMS error.

clear all; close all; beep off; clc;                                        % Close and clear all

%% Input needed for MPM
constant = struct('density',1e3,'E',1e5,'g',-9.81,'load',0,...              % ...
                  'height',1,'alpha',0);                                    % Define constants
              
% Define settings of the flags (MPM, FEM, ULFEM)
% MPM   = [0 1 1 1 1 0 1 1 0 0 1 0 0 0]
% FEM   = [0 0 1 0 0 1 0 0 0 0 0 0 0 0]
% ULFEM = [0 0 1 1 0 1 0 1 0 0 0 1 0 0]
flag     = struct('both_ends_fixed',0,'volume_update',1,'lumped',1,'change_glob_pos',1,'change_loc_pos',1,...
                'lagranian',0,'momentum',1,'deformation',1,'num_int',0,'dynamic',0,'splines',1,'ULFEM',0);
       
%% Knot vector and it's properties
deg = 2;                                                                    % Degree of basis functions    
n_dof = 66;                                                                 % Number of degrees of freedom (DOF)
Xi = zeros(n_dof + deg + 1,1);                                              % Knot vector initialized
mesh = 0:constant.height/(n_dof-deg):constant.height;                       % ...
Xi(deg+2:n_dof+1,1) = mesh(2:n_dof-deg+1);                                  % ...
Xi(end-deg:end,1) = mesh(end);                                              % Knot vector determined
knot_partition = Xi(2:end,1) - Xi(1:end-1);                                 % Distance between knots
number_knot_elements = nnz(knot_partition);                                 % Number of (non-zero) knot spans
min_knot_element_size = min(nonzeros(Xi));                                  % Minimal knot span size (non-zero)

%% Particle properties
number_particles_per_knot_element = 4;                                      % Number of particles per knot span
number_particles = number_knot_elements*number_particles_per_knot_element;  % Total number of particles 

if flag.num_int == 1
    loc = [1/4 - (1/4)*(1/sqrt(3)),1/4 + (1/4)*(1/sqrt(3)), ...             % ...
           3/4 - (1/4)*(1/sqrt(3)), 3/4 + (1/4)*(1/sqrt(3))];               % Local position integration points/particles       
    weight = constant.height/number_particles*ones...                       % ...
            (number_knot_elements*number_particles_per_knot_element,1);     % Weight integration points/particles
else
    loc = [1:number_particles_per_knot_element]/...                         % ...
          (number_particles_per_knot_element+1);                            % ...
    loc = [1:2:(2*number_particles_per_knot_element-1)]/...                 % ...
          (2*number_particles_per_knot_element);                            % Local position integration point/particles
    weight = min_knot_element_size/number_particles_per_knot_element*...    % ...
        ones(number_knot_elements*number_particles_per_knot_element,1);     % Weight integration points/particles
end

pos_p_loc  = repmat(loc,1,number_knot_elements);                            % Local position particles  
pos_p_glob = pos_p_loc'.*kron(nonzeros(knot_partition),ones...              % ...
             (number_particles_per_knot_element,1)) + kron(Xi(1+deg:...     % ...
             end-deg-1), ones(number_particles_per_knot_element,1));        % Global position particles
pos_p_glob(end) = 0.995;
%% Time step size, etc.
CFL_number = 0.95;                                                          % CFL number
total_time = 2.5;                                                           % Total time simulation
t_cr = min_knot_element_size/sqrt(constant.E/constant.density);             % Critical time step
t_step = 1e-4;                                                              % Time step size
number_time_steps = floor(total_time/t_step);                               % Number of time steps 
t = 0:t_step:(number_time_steps-1)*t_step;                                  % Time vector

%% Initial conditions 
displacement_func = @(x) 0*x;                                               % Initial displacement
velocity_func = @(x) 0*x;                                                   % Initial velocity
stress_func = @(x) 0*x;                                                     % Initial stress

%% Energy vectors
E_kin = zeros(1,number_time_steps);                                         % Kinetic energy
E_pot = zeros(1,number_time_steps);                                         % Potential energy
E_grav = zeros(1,number_time_steps);                                        % Gravitational energy
E_trac = zeros(1,number_time_steps);                                        % Traction energy
Total_F = zeros(n_dof,number_time_steps);                                   % Total force

%% Compute the solution using MPM
tic
[displacement_mpm, velocity_mpm,velocity_mpm_nodes, M_lump,...
    displacement_mpm_particles, E_kin,E_pot,E_grav,E_trac, ...
    stress_p, strain_p] = MPM_1D_B_spline(constant,flag,pos_p_glob, pos_p_loc, number_knot_elements, min_knot_element_size,...
    number_particles_per_knot_element, t_step, number_time_steps, total_time, ...
    E_kin,E_pot,E_grav, E_trac,weight,Xi,deg,displacement_func,velocity_func,stress_func,number_particles,Total_F);
toc
position_mpm_particles = displacement_mpm_particles + repmat(pos_p_glob,1, number_time_steps); % Position particles

%% Obtain the exact solution for the particles(!)
for p = 1:number_particles
    [position_exact(p,:) sol_exact(p,:) vel_exact(p,:)] = exact_solution...
        (constant.density,constant.E,constant.load,-constant.g,...
        constant.height,pos_p_glob(p), t);
end

%% Plot displacement/velocity versus time for one particle
particle_number = floor(1*number_particles/2);                              % particle number

figure(1)
set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);
plot(t,position_mpm_particles(particle_number,:), 'LineWidth',2)
hold on
plot(t,position_exact(particle_number,:),'--r', 'LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('position [m]','FontSize', 12)
title(sprintf('Position of particle %d',particle_number),'FontSize', 12)
if flag.change_loc_pos == 0
    legend('FEM','Exact')
else
    legend('MPM','Exact')
end
hold on

figure(2)
set(gcf, 'PaperPosition', [0 0 7 7]);
set(gcf, 'PaperSize', [6 4.5]);
plot(t,velocity_mpm(particle_number,:),'b','LineWidth',2)
hold on
plot(t,vel_exact(particle_number,:),'--r','LineWidth',2)
xlabel('time [s]','Fontsize',12)
ylabel('velocity [m/s]', 'Fontsize',12)
title(sprintf('Velocity of particle %d',particle_number),'FontSize', 12)
if flag.change_loc_pos == 0
    legend('FEM','Exact')
else
    legend(' MPM','Exact')
end

E_total = - E_kin + E_pot + E_trac + E_grav;
figure(3)
set(gcf, 'PaperPosition', [0 0 7 7]);
set(gcf, 'PaperSize', [6 4.5]);
plot(t,-E_kin,'-b','LineWidth',2)
hold on
plot(t,E_pot,'-r','LineWidth',2)
plot(t,E_grav,'-g','LineWidth',2)
plot(t,E_trac,'-k','LineWidth',2)
plot(t,E_total,'-m','LineWidth',2)
xlabel('Time [s]','Fontsize',12)
ylabel('Energy','Fontsize',12)
legend('E_{kin}','E_{pot}' ,'E_{grav}', 'E_{trac}', 'E_{total}')

%% Accuracy of MPM solution
t_check = floor(length(t)/5);                                               % Time of comparison 
Error = norm(position_exact(:,floor(t_check))-position_mpm_particles(:,...  % ...
    floor(t_check)))/sqrt(number_particles)                                 % RMS error




