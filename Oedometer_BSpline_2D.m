%% Oedometer_Bspline_2D (Bottom end fixed) 
%Pascal de Koster, based on code by Roel Tielen and Lisa Wobbes
%This code aims to solve the motion of a vibrating bar with both ends
%fixed, where the bar is simulated as a two dimensional object, although
%this problem can be posed as one dimensional. The purpose is to validate
%the 2 dimensional code for this relatively easy scenario. Note that the
%height in this code is the x-direction, so x and y are swapped compared to
%the report. Only in the final figures y is the height.

clear; close all;% clc;
profile on
tic

%% Input constants

%The constants of the problem: from left to right:
%Dimension, density, Young's modulus, gravitational constant, the presence
%of a load (external gravitational force) on the top of the structure, the
%length
% constant = struct('dim',2,'density',1,'E',5e4,'g',-9.81e0,'load',0,...
%                   'length',25,'width',2,'alpha',0,'v_0',0,...
%                   'Poisson_ratio',0);
constant = struct('dim',2,'density',1e3,'E',1e5,'g',-9.81e0,'load',0,...
                  'length',1,'width',0.1,'alpha',0,'v_0',0,...
                  'Poisson_ratio',0);

flag = struct('spline_basis',1,'both_ends_fixed',0,'lumped',2,...
    'Gauss_integration',0,'TLS',1,'Cubic_TLS',0,'grid_regularity',1,...
    'symmetric',1,'particle_regularity',1,'DisableNodes',0);

%Set the seed for random number to obtain the same result every time
rng(2)

%% Triangular grid and its properties

%Number of vertices in x- and y-direction
n_vert_X = 16+1;  
n_vert_Y = 2+1;

%Number of ghost verticex in x- and y-direction
n_vert_X_ghost = 1;
n_vert_Y_ghost = 0;

%Total number of vertices
n_vertices_X = n_vert_X + n_vert_X_ghost;
n_vertices_Y = n_vert_Y + n_vert_Y_ghost;

%The regularity of the grid, 0 is for random, 1 for regular.
grid_regularity=flag.grid_regularity;
%Symmetric only if grid_regularity = 1
symmetric=flag.symmetric;
if grid_regularity == 0
    symmetric = 0;
end

if symmetric==1
    n_vertices=n_vertices_X*n_vertices_Y+(n_vertices_X-1)*(n_vertices_Y-1);
else
    n_vertices=n_vertices_X*n_vertices_Y;
end

%number of degrees of freedom
if flag.spline_basis==1
    n_dof = (constant.dim+1)*n_vertices;
else
    n_dof = n_vertices;
end

%Length with ghost points added
ghost_length = constant.length*(n_vertices_X-1)/(n_vert_X-1);
    
if grid_regularity==1
    %Vertices in rectangular order
    vertices_X=linspace(0,ghost_length,n_vertices_X)'...
        *ones(1,n_vertices_Y); 
    vertices_X=reshape(vertices_X,[n_vertices_X*n_vertices_Y,1]);
    if symmetric ==1
        %If symmetric, then add center nodes in each rectangle to make the
        %grid symmetric.
        vertices_X2=ghost_length*linspace(0.5/(n_vertices_X-1),...
            1-0.5/(n_vertices_X-1),n_vertices_X-1)'*...
            ones(1,n_vertices_Y-1);
        vertices_X2=reshape(vertices_X2,...
            [(n_vertices_X-1)*(n_vertices_Y-1),1]);
        vertices_X=[vertices_X;vertices_X2];
    end
    vertices_Y=ones(n_vertices_X,1)...
        *linspace(0,constant.width,n_vertices_Y);
    vertices_Y=reshape(vertices_Y,[n_vertices_X*n_vertices_Y,1]);
    if symmetric ==1
        %If symmetric, then add center nodes in each rectangle to make the
        %grid symmetric.
        vertices_Y2=ones(n_vertices_X-1,1)*...
            constant.width*linspace(0.5/(n_vertices_Y-1),...
            1-0.5/(n_vertices_Y-1),n_vertices_Y-1);            
        vertices_Y2=reshape(vertices_Y2,...
            [(n_vertices_X-1)*(n_vertices_Y-1),1]);
        vertices_Y=[vertices_Y;vertices_Y2];
    end
    
    %Vertex Numbering goes as follows:
    %13 14 15 16
    % 9 10 11 12
    % 5  6  7  8
    % 1  2  3  4
    
elseif grid_regularity == 2     %ordened with random displacements
    rand_dis_max_frac = 0.2; %The fraction of extra displacement wrt 
                             %the element distance (<0.5)
    %The maximum added random displacement vertices in rectangular order
    rand_dis_max_X = rand_dis_max_frac * constant.length/(n_vert_X-1); 
    rand_dis_max_Y = rand_dis_max_frac * constant.width/(n_vert_Y-1); 
    
    vertices_X=linspace(0,ghost_length,n_vertices_X)'...
        *ones(1,n_vertices_Y);
    %Add randomness to internal particles
    vertices_X = vertices_X +...
        [zeros(1,n_vertices_Y);...
         (rand(n_vertices_X-2,n_vertices_Y)-0.5)*rand_dis_max_X;...
         zeros(1,n_vertices_Y)];
     
    vertices_X=reshape(vertices_X,[n_vertices_X*n_vertices_Y,1]);
    if symmetric == 1
        %If symmetric, then add center nodes in each rectangle to make the
        %grid symmetric.
        vertices_X2=ghost_length*linspace(0.5/(n_vertices_X-1),...
            1-0.5/(n_vertices_X-1),n_vertices_X-1)'*...
            ones(1,n_vertices_Y-1)+...
            (rand(n_vertices_X-1,n_vertices_Y-1)-0.5)*rand_dis_max_X;
        vertices_X2=reshape(vertices_X2,...
            [(n_vertices_X-1)*(n_vertices_Y-1),1]);
        vertices_X=[vertices_X;vertices_X2];
    end
    vertices_Y=ones(n_vertices_X,1)...
        *linspace(0,constant.width,n_vertices_Y);
    %Add randomness to internal particles
    vertices_Y = vertices_Y +...
        [zeros(n_vertices_X,1),...
         (rand(n_vertices_X,n_vertices_Y-2)-0.5)*rand_dis_max_Y,...
         zeros(n_vertices_X,1)];
     
    vertices_Y=reshape(vertices_Y,[n_vertices_X*n_vertices_Y,1]);
    if symmetric ==1
        %If symmetric, then add center nodes in each rectangle to make the
        %grid symmetric.
        vertices_Y2=ones(n_vertices_X-1,1)*...
            constant.width*linspace(0.5/(n_vertices_Y-1),...
            1-0.5/(n_vertices_Y-1),n_vertices_Y-1)+...
            (rand(n_vertices_X-1,n_vertices_Y-1)-0.5)*rand_dis_max_Y;            
        vertices_Y2=reshape(vertices_Y2,...
            [(n_vertices_X-1)*(n_vertices_Y-1),1]);
        vertices_Y=[vertices_Y;vertices_Y2];
    end
else
    %Irregular grid
    %The center part is randomly generated, the edges are regular
    vertices_X=[ zeros(1,n_vertices_Y-2);...
                 ghost_lengthh*rand(n_vertices_X-2,n_vertices_Y-2);...
                 ghost_length*ones(1,n_vertices_Y-2) ];
    vertices_X=[ linspace(0,ghost_length,n_vertices_X)',...
                 vertices_X,...
                 linspace(0,ghost_length,n_vertices_X)'];
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
trimesh(triangles,vertices_Y,vertices_X)
hold off
% title('Triangulations')

%% Time step size, etc.
CFL_number = 0.95;                                                          % CFL number
total_time = 2.5 ;                                                          % Total time simulation
% total_time = 0.5;                                                       % Total time simulation
% t_cr = min_knot_element_size/sqrt(constant.E/constant.density);           % Critical time step
t_step = 1e-3;                                                            % Time step size
n_time_steps = floor(total_time/t_step)+1;                                  % Number of time steps 
t = 0:t_step:(n_time_steps-1)*t_step;                                       % Time vector

%% Particle properties

%The regularity of the particles, 0 is for completely random, 1 for
%completely regular.
particle_regularity=flag.particle_regularity;

%Number of particles for each direction on the boundaries
%IMPORTANT NOTE! In case of regular particles and regular grid, it often
%occurs that particles lie at the start exactly on a boundary between two
%triangulations. Try to prevent this in order to avoid unnecessary
%inaccuracies.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Number of x- and y-particle EVEN and DIFFERENT, to prevent difficulties
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n_particles_per_node_X = 14;                                            
n_particles_per_node_Y = 12;
n_particles_X = n_particles_per_node_X * (n_vert_X-1);
n_particles_Y = n_particles_per_node_Y * (n_vert_Y-1);
%Number of particles spreaded over the grid
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
elseif particle_regularity==2
    rand_dis_max_frac_p = 0.1;
    
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
    
    particles_X=particles_X+rand_dis_max_frac_p*...
        constant.length/n_particles_X*(rand(size(particles_X))-0.5);
    particles_Y=particles_Y+rand_dis_max_frac_p*...
        constant.width/n_particles_Y*(rand(size(particles_Y))-0.5);
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
scatter(particles_Y(:,1),particles_X(:,1),10,'r')
% title('Triangulations and material points')


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
% scatter(particles_X(:,1),particles_Y(:,1),'r')
% title('Material points and weigth surfaces')

%% Initial conditions 

%These conditions are identical conditions of the 1D case, extended in the
%y-direction.
displacement_func = @(x,y) 0*x;                                             % Initial displacement

% Constant for exact solution
w1 = pi*sqrt(constant.E/constant.density)/(constant.length);                % Constant for exact solution
b1 = pi/(constant.length);   

velocity_Y_func_exact = @(x,y,t) 0.*x;                                      % Initial velocity Y
stress_func = @(x,y) 0.*x;                                                  % Initial stress

velocity_X_func = @(x,y) 0.*x.*y;
velocity_Y_func = @(x,y) 0.*x.*y;

%% Obtain the exact solution for the particles(!)

position_X_exact=zeros(n_particles_X,n_time_steps);
disp_X_exact = zeros(n_particles_X,n_time_steps);
disp_Y_exact = zeros(n_particles,n_time_steps);
vel_X_exact  = zeros(n_particles_X,n_time_steps);
vel_Y_exact  = zeros(n_particles,n_time_steps);

for p = 1:n_particles_X
    [position_X_exact(p,:), disp_X_exact(p,:), vel_X_exact(p,:)] = ...
        exact_solution_Oedometer...
        (constant.density,constant.E,constant.load,-constant.g,...
        constant.length,particles_X(p), t);
end
    
position_X_exact = repmat(position_X_exact,[n_particles_Y,1]);
disp_X_exact = repmat(disp_X_exact,[n_particles_Y,1]);
vel_X_exact = repmat(vel_X_exact,[n_particles_Y,1]);

%% Compute the solution using MPM

[u_X_n, u_Y_n, particles_X, particles_Y, v_X_p, v_Y_p,...
    v_X_n, v_Y_n, mass_p,u_X_p, u_Y_p,...
    E_kin,E_pot,E_grav,E_trac,...
    stress11_p,stress12_p,stress21_p,stress22_p, ...
    strain11_p,strain12_p,strain21_p,strain22_p, ...
    N_vec, B_vec_X, B_vec_Y, volume_p] = ...
    MPM_2D_B_spline(constant,flag,n_dof,n_triangles,...
    n_particles,particles_X,particles_Y,vertices_X,vertices_Y,triangles,...
    area_weight,velocity_X_func,velocity_Y_func,t_step,...
    n_time_steps);


%% Plotting results

%Choose exact = 1 for small deformations, and exact = 0 for large
%deformations and use the ULFEM solution
exact = 0;
if exact == 0
    load('Results\SC_LD_ULFEM.mat')
end

prt_center = floor(1/2*n_particles_Y)*n_particles_X+...
             floor(1/2*n_particles_X);

%Displacement X
figure
hold on
plot(t,particles_X(prt_center,:)-...
    particles_X(prt_center,1),...
    'b','LineWidth',2)
if exact == 1
    plot(t,position_X_exact(prt_center,:)-...
        particles_X(prt_center,1),...
        'r--','LineWidth',2)
    legend('MPM','exact')
else %ULFEM
    plot(ULFEM_t,ULFEM_disp_center,...
        'r--','LineWidth',2)
    legend('MPM','ULFEM')
end
% title('X-displacement')
xlabel('t [s]')
ylabel('Y-displacement of middle particle [m]')
hold off

%Velocity Y
figure
hold on
plot(t,v_Y_p(prt_center,:),...
    'b','LineWidth',2)
plot(t,vel_Y_exact(prt_center,:),...
    'r--','LineWidth',2)
% title('Y-velocity')
xlabel('t [s]')
ylabel('X-velocity of middle particle [m/s]')
legend('MPM','exact')
hold off

%Velocity X
figure
hold on
plot(t,v_X_p(prt_center,:),...
    'b','LineWidth',2)
if exact == 1
    plot(t,vel_X_exact(prt_center,:),...
        'r--','LineWidth',2)
    legend('MPM','ULFEM')
else %ULFEM
    plot(ULFEM_t,ULFEM_velocity_center,...
        'r--','LineWidth',2)
    legend('MPM','ULFEM')
end
% title('X-velocity')
xlabel('t [s]')
ylabel('Y-velocity of middle particle [m/s]')

hold off

figure
hold on
plot(stress11_p(1:n_particles_X,end),particles_X(1:n_particles_X,1),...
    'b','LineWidth',2)
if exact == 1
else %ULFEM
    plot(ULFEM_t,ULFEM_velocity_center,...
        'r--','LineWidth',2)
end
% title('X-velocity')
xlabel('Initial height [m]')
ylabel('normal stress of middle particle in y-direction [m/s]')
% legend('MPM','exact')
hold off


fprintf('Dt=%.1e; (n_vx,nv_y)=(%i,%i); (n_px,n_py)=(%i,%i)  \n',...
    t_step,n_vertices_X,n_vertices_Y,...
    n_particles_per_node_X,n_particles_per_node_Y);

fprintf('The error at t_end for the middle particle is %.1e \n',...
    v_X_p(floor((n_vertices_X-1)*n_particles_per_node_X/2),end)-...
    vel_X_exact(floor((n_vertices_X-1)*n_particles_per_node_X/2),end) );

fprintf('The error integrated over time up to t_end is %.1e \n',...
    t_step*sum(abs( ...
    v_X_p(floor((n_vertices_X-1)*n_particles_per_node_X/2),:)-...
    vel_X_exact(floor((n_vertices_X-1)*n_particles_per_node_X/2),:))) );

%% L2-Error at t_end

figure
hold on
plot(particles_X(1:n_particles_X,1),...
    v_X_p(1:n_particles_X,n_time_steps),'b')
plot(particles_X(1:n_particles_X,1),...
    vel_X_exact(1:n_particles_X,n_time_steps),'r--')
title(sprintf('Velocity at t=%.1e',total_time))
ylabel('v [m/s]')
xlabel('y [m]')
legend('MPM','exact')

figure
plot(particles_X(1:n_particles_X,1),...
    (v_X_p(1:n_particles_X,n_time_steps)-...
    vel_X_exact(1:n_particles_X,n_time_steps)) )
title('Displacement error')
ylabel('error [m]')
xlabel('y [m]')

E_vel=sqrt(sum(volume_p(:,n_time_steps).*(v_X_p(:,n_time_steps)-...
    vel_X_exact(:,n_time_steps)).^2));
fprintf('E_vel = %0.5e \n', E_vel);

toc
profile viewer


