%% Vib_bar_Bspline_2D (both ends fixed) 
%Pascal de Koster, based on code by Roel Tielen and Lisa Wobbes
%This code aims to solve the motion of a vibrating bar with both ends
%fixed, where the bar is simulated as a two dimensional object, although
%this problem can be posed as one dimensional. The purpose is to validate
%the 2 dimensional code for this relatively easy scenario.

clear; close all; %clc;
% profile on
%
%% Input constants

%The constants of the problem: from left to right:
%Dimension, density, Young's modulus, gravitational constant, the presence
%of a load (external gravitational force) on the top of the structure, the
%length
constant = struct('dim',2,'density',1,'E',100,'g',0,'load',0,...
                  'length',25,'width',5,'alpha',0,'v_0',0.1,...
                  'Poisson_ratio',0);
    
flag     = struct('spline_basis',1,'both_ends_fixed',1,'lumped',0,'splines',0);

%% Triangular grid and its properties

deg = 2;                                                                    % Degree of basis functions  
n_vertices_X = 32;
n_vertices_Y = 4;
n_vertices = n_vertices_X*n_vertices_Y;                                     % Number of degrees of freedom (DOF)

if flag.spline_basis==1
    n_dof = (constant.dim+1)*n_vertices;
else
    n_dof = n_vertices;
end

%The regularity of the grid, 0 is for random, 1 for regular.
grid_regularity=1;

if grid_regularity==1
    %Vertices in rectangular order 
    vertices_X=linspace(0,constant.length,n_vertices_X)'...
        *ones(1,n_vertices_Y);
    vertices_X=reshape(vertices_X,[n_vertices_X*n_vertices_Y,1]);
    vertices_Y=ones(n_vertices_X,1)...
        *linspace(0,constant.width,n_vertices_Y);
    vertices_Y=reshape(vertices_Y,[n_vertices_X*n_vertices_Y,1]);
    %Vertex Numbering goes as follows:
    %13 14 15 16
    % 9 10 11 12
    % 5  6  7  8
    % 1  2  3  4
else
    %Irregular grid
    %The center part is randomly generated, the edges are regular
    vertices_X=[ zeros(1,n_vertices_Y-2);...
                 constant.length*rand(n_vertices_X-2,n_vertices_Y-2);...
                 constant.length*ones(1,n_vertices_Y-2) ];
    vertices_X=[ linspace(0,constant.length,n_vertices_X)',...
                 vertices_X,...
                 linspace(0,constant.length,n_vertices_X)'];
    vertices_X=reshape(vertices_X,[n_vertices_X*n_vertices_Y,1]);
    
    vertices_Y=[ zeros(n_vertices_X-2,1),...
                 constant.width*rand(n_vertices_X-2,n_vertices_Y-2),...
                 constant.width*ones(n_vertices_X-2,1)];
    vertices_Y=[ linspace(0,constant.width,n_vertices_Y);...
                 vertices_Y;...
                 linspace(0,constant.width,n_vertices_Y) ];
    vertices_Y=reshape(vertices_Y,[n_vertices_X*n_vertices_Y,1]);
end

%Generate triangles
triangles = delaunayTriangulation(vertices_X,vertices_Y);
n_triangles=size(triangles,1);

figure(1)
hold on
trimesh(triangles,vertices_X,vertices_Y)
hold off
title('Triangulations')

%% Time step size, etc.
CFL_number = 0.95;                                                          % CFL number
total_time = 2.5;                                                           % Total time simulation
% total_time = 2e-2;                                                        % Total time simulation
% t_cr = min_knot_element_size/sqrt(constant.E/constant.density);           % Critical time step
t_step = 1e-2;                                                              % Time step size
n_time_steps = floor(total_time/t_step);                                    % Number of time steps 
t = 0:t_step:(n_time_steps-1)*t_step;                                       % Time vector

%% Particle properties

%The regularity of the particles, 0 is for completely random, 1 for
%completely regular.
particle_regularity=1;

%Number of particles for each direction on the boundaries
%IMPORTANT NOTE! In case of regular particles and regular grid, it often
%occurs that particles lie at the start exactly on a boundary between two
%triangulations. Try to prevent this in order to avoid unnecessary
%inaccuracies.
n_particles_per_node_X = 5;                                            
n_particles_per_node_Y = 4;
n_particles_X = n_particles_per_node_X * (n_vertices_X-1);
n_particles_Y = n_particles_per_node_Y * (n_vertices_Y-1);
%Number of particles spreaded of the grid
n_particles   = n_particles_X * n_particles_Y;

particles_X=zeros(n_particles,n_time_steps);
particles_Y=zeros(n_particles,n_time_steps);

if particle_regularity==1
    particles_X_temp=(1/2:n_particles_X-1/2)'/...
        n_particles_X*constant.length;
    particles_X_temp=particles_X_temp*ones(1,n_particles_Y);
    particles_X(:,1)=reshape(particles_X_temp,...
        [n_particles_X*n_particles_Y,1]);
    
    particles_Y_temp=(1/2:n_particles_Y-1/2)/...
        n_particles_Y*constant.width;
    particles_Y_temp=ones(n_particles_X,1)*particles_Y_temp;
    particles_Y(:,1)=reshape(particles_Y_temp,...
        [n_particles_X*n_particles_Y,1]);
else
    particles_X(:,1)=constant.length*rand(n_particles,1);
    particles_Y(:,1)=constant.width *rand(n_particles,1);
end

%Properties, the area of each particle, inlcuding the weight
[area_vertex_X,area_vertex_Y,area_vertex_indices,area_weight] = ...
    material_point_area(constant,particles_X(:,1),particles_Y(:,1));  

%plot the material points and the corresponding surfaces
figure(1)
hold on
scatter(particles_X(:,1),particles_Y(:,1))
title('Triangulations and material points')


%Uncomment for illustrating the weight of each particle 
%
% figure
% hold on
% for i=1:length(area_vertex_indices)
%     ind_temp=[cell2mat(area_vertex_indices(i)),0];
%     ind_temp(end)=ind_temp(1);
%     plot(area_vertex_X(ind_temp),area_vertex_Y(ind_temp),'LineWidth',...
%         3*area_weight(i))
% end
% title('Material points and weigth surfaces')

%% Initial conditions 

%These conditions are identical conditions of the 1D case, extended in the
%y-direction.
displacement_func = @(x,y) 0*x;                                             % Initial displacement

% Constant for exact solution
w1 = pi*sqrt(constant.E/constant.density)/(constant.length);                % Constant for exact solution
b1 = pi/(constant.length);   

velocity_X_func_exact = @(x,y,t) constant.v_0*...
    sin(pi*x/((1+(1-flag.both_ends_fixed))*constant.length))*cos(w1*t);     % Exact velocity X
vel_X_exact(:,1) = velocity_X_func_exact(particles_X(:,1),...
        particles_Y(:,1),0);                                                % Initial velocity X
velocity_Y_func_exact = @(x,y,t) 0.*x;                                      % Initial velocity Y
stress_func = @(x,y) 0.*x;                                                  % Initial stress

%% Obtain the exact solution for the particles(!)

disp_X_exact = zeros(n_particles,n_time_steps);
vel_Y_exact  = zeros(n_particles,n_time_steps);
disp_Y_exact = zeros(n_particles,n_time_steps);

for n = 1:n_time_steps-1
    vel_X_exact(:,n+1) = velocity_X_func_exact(particles_X(:,1),...
        particles_Y(:,1),n*t_step);
    disp_X_exact(:,n+1) = constant.v_0/w1*sin(w1*t_step*n)*...
        sin(b1*particles_X(:,1));
end
position_exact = disp_X_exact + repmat(particles_X(:,1),1,n_time_steps);    % Exact position particles

%% Compute the solution using MPM
tic
[u_X_n, u_Y_n, particles_X, particles_Y, v_X_p, v_Y_p,...
    v_X_n, v_Y_n, mass_p,u_X_p, u_Y_p,...
    E_kin,E_pot,E_grav,E_trac,...
    stress_p,strain_p,N_vec, B_vec_X, B_vec_Y] = ...
    MPM_2D_B_spline(constant,flag,n_dof,n_triangles,...
    n_particles,particles_X,particles_Y,vertices_X,vertices_Y,triangles,...
    area_weight,velocity_X_func_exact,velocity_Y_func_exact,t_step,...
    n_time_steps);
toc

%% Plotting results

%Displacement X
figure
hold on
plot(t,particles_X(floor(n_vertices_X*n_particles_per_node_X/2),:),'b')
plot(t,position_exact(floor(n_vertices_X*n_particles_per_node_X/2),:),'r')
title('X-displacement')
xlabel('t')
ylabel('Position of middle particle')
legend('MPM','exact')
hold off

%Velocity Y
figure
hold on
plot(t,v_Y_p(floor(n_vertices_X*n_particles_per_node_X/2),:),'b')
plot(t,vel_Y_exact(floor(n_vertices_X*n_particles_per_node_X/2),:),'r')
title('Y-velocity')
xlabel('t')
ylabel('Position of middle particle')
legend('MPM','exact')
hold off

%Velocity X
figure
hold on
plot(t,v_X_p(floor(n_vertices_X*n_particles_per_node_X/2),:),'b')
plot(t,vel_X_exact(floor(n_vertices_X*n_particles_per_node_X/2),:),'r')
title('X-velocity')
xlabel('t')
ylabel('Position of middle particle')
legend('MPM','exact')
hold off

fprintf('Dt=%.1e; (n_vx,nv_y)=(%i,%i); (n_px,n_py)=(%i,%i)  \n',...
    t_step,n_vertices_X,n_vertices_Y,...
    n_particles_per_node_X,n_particles_per_node_Y);

fprintf('The error at t_end for the middle particle is %.1e \n',...
    v_X_p(floor(n_vertices_X*n_particles_per_node_X/2),end)-...
    vel_X_exact(floor(n_vertices_X*n_particles_per_node_X/2),end) );

fprintf('The error integrated over time up to t_end is \n')
fprintf('    %.1e \n', t_step*sum(abs( ...
    v_X_p(floor(n_vertices_X*n_particles_per_node_X/2),:)-...
    vel_X_exact(floor(n_vertices_X*n_particles_per_node_X/2),:))) );

% profile viewer


