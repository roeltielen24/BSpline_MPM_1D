%% Vib_bar_Bspline_2D (both ends fixed) 
%Pascal de Koster, based on code by Roel Tielen and Lisa Wobbes
%This code aims to solve the motion of a vibrating bar with both ends
%fixed, where the bar is simulated as a two dimensional object, although
%this problem can be posed as one dimensional. The purpose is to validate
%the 2 dimensional code for this relatively easy scenario.

clear; close all;% clc;
profile on
tic

%% Input constants

%The constants of the problem: from left to right:
%Dimension, density, Young's modulus, gravitational constant, the presence
%of a load (external gravitational force) on the top of the structure, the
%length
%Small deformations
constant = struct('dim',2,'density',1,'E',100,'g',0,'load',0,...
                  'length',25,'width',2,'alpha',0,'v_0',0.1,...
                  'Poisson_ratio',0);
% %Large deformations
% constant = struct('dim',2,'density',25,'E',50,'g',0,'load',0,...
%                   'length',1,'width',0.1,'alpha',0,'v_0',0.1,...
%                   'Poisson_ratio',0);

flag = struct('spline_basis',0,'both_ends_fixed',1,'lumped',1,...
    'Gauss_integration',0,'TLS',1,'Cubic_TLS',0,'grid_regularity',1,...
    'symmetric',1,'particle_regularity',1,'DisableNodes',0);

%Set the seed for random number to obtain the same result every time
rng(2)

%% Triangular grid and its properties

% Degree of basis functions 
deg = 2;
n_vertices_X = 24+1;
n_vertices_Y = 1+1;

%The regularity of the grid, 0 is for random, 1 for regular.
grid_regularity=flag.grid_regularity;
%Symmetric only if grid_regularity = 1
symmetric=flag.symmetric;
if grid_regularity ~= 1
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

if grid_regularity==1
    %Vertices in rectangular order 
    vertices_X=linspace(0,constant.length,n_vertices_X)'...
        *ones(1,n_vertices_Y); 
    vertices_X=reshape(vertices_X,[n_vertices_X*n_vertices_Y,1]);
    if symmetric ==1
        %If symmetric, then add center nodes in each rectangle to make the
        %grid symmetric.
        vertices_X2=constant.length*linspace(0.5/(n_vertices_X-1),...
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
% title('Triangulations')

%% Time step size, etc.
CFL_number = 0.95;                                                          % CFL number
total_time = 2.5 ;                                                           % Total time simulation
% total_time = 2e-2;                                                        % Total time simulation
% t_cr = min_knot_element_size/sqrt(constant.E/constant.density);           % Critical time step
t_step = 1e-3;                                                              % Time step size
n_time_steps = floor(total_time/t_step)+1;                                    % Number of time steps 
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
%Number of x- and y-particle should be EVEN and DIFFERENT, to prevent
%difficulties Suggested: n_x = 16, n_y = 14.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n_particles_per_node_X = 10;                                            
n_particles_per_node_Y = 8;
n_particles_X = n_particles_per_node_X * (n_vertices_X-1);
n_particles_Y = n_particles_per_node_Y * (n_vertices_Y-1);
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
else
    particles_X(:,1)=constant.length*rand(n_particles,1);
    particles_Y(:,1)=constant.width *rand(n_particles,1);
end

%Properties, the area of each particle, inlcuding the weight
[area_vertex_X,area_vertex_Y,area_vertex_indices,area_weight] = ...
    material_point_area(constant,particles_X(:,1),particles_Y(:,1));  

% %Uncomment for plotting the material points and the corresponding surfaces
figure(1)
hold on
scatter(particles_X(:,1),particles_Y(:,1),3,'r')
% title('Triangulations and material points')


% %Uncomment for illustrating the weight of each particle 
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

vel_X_exact  = zeros(n_particles,n_time_steps);

velocity_X_func_exact = @(x,y,t) constant.v_0*...
    sin(pi*x/((1+(1-flag.both_ends_fixed))*constant.length))*cos(w1*t);     % Exact velocity X
vel_X_exact(:,1) = velocity_X_func_exact(particles_X(:,1),...
        particles_Y(:,1),0);                                                % Initial velocity X
velocity_Y_func_exact = @(x,y,t) 0.*x;                                      % Initial velocity Y
stress_func = @(x,y) 0.*x;                                                  % Initial stress

velocity_X_func = @(x,y) constant.v_0*...
    sin(pi*x/((1+(1-flag.both_ends_fixed))*constant.length));
velocity_Y_func = @(x,y) 0*x.*y;

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
exact = 1;
if exact == 0
    load('Results\VB_LD_ULFEM.mat')
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
    plot(t,position_exact(prt_center,:)-...
        particles_X(prt_center,1),...
        'r--','LineWidth',2)
else %ULFEM
    plot(ULFEM_t,ULFEM_disp_center,...
        'r--','LineWidth',2)
end
% title('X-displacement')
xlabel('time [s]')
ylabel('X-displacement of middle particle [m]')
if exact == 1
    legend('MPM','exact')
else %ULFEM
    legend('MPM','ULFEM')
end
hold off

%Velocity Y
figure
hold on
plot(t,v_Y_p(prt_center,:),...
    'b','LineWidth',2)
plot(t,vel_Y_exact(prt_center,:),...
    'r--','LineWidth',2)
% title('Y-velocity [m/s]')
xlabel('time [s]')
ylabel('Y-velocity of middle particle [m/s]')
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
else %ULFEM
    plot(ULFEM_t,ULFEM_velocity_center,...
        'r--','LineWidth',2)
end
% title('X-velocity')
xlabel('time [s]')
ylabel('X-velocity of middle particle [m/s]')
if exact == 1
    legend('MPM','exact')
else %ULFEM
    legend('MPM','ULFEM')
end
hold off

fprintf('Dt=%.1e; (n_vx,nv_y)=(%i,%i); (n_px,n_py)=(%i,%i)  \n',...
    t_step,n_vertices_X,n_vertices_Y,...
    n_particles_per_node_X,n_particles_per_node_Y);

fprintf('The error at t_end for the middle particle is %.1e \n',...
    v_X_p(prt_center,end)-...
    vel_X_exact(prt_center,end) );

fprintf('The error integrated over time up to t_end is %.1e \n',...
    t_step*sum(abs( ...
    v_X_p(prt_center,:)-...
    vel_X_exact(prt_center,:))) );

%% L2-Error at t_end

%Stress X
figure
hold on
plot(particles_X(1:n_particles_X,1),...
     stress11_p(1:n_particles_X,n_time_steps))
% title(sprintf('\\sigma_{11} at t=%.1e',total_time))
xlabel('x')
ylabel('\sigma_{11}')
hold off


figure
hold on
plot(particles_X(1:n_particles_X,1),...
    v_X_p(1:n_particles_X,n_time_steps),'b')
plot(particles_X(1:n_particles_X,1),...
    velocity_X_func_exact(particles_X(1:n_particles_X,1),...
    particles_Y(1:n_particles_X,1),total_time),'r--')
% title(sprintf('Velocity at t=%.1e',total_time))
ylabel('velocity [m/s]')
xlabel('x')
legend('MPM','exact')


figure
plot(particles_X(1:n_particles_X,1),...
    (u_X_p(1:n_particles_X,n_time_steps)-...
    disp_X_exact(1:n_particles_X,n_time_steps)) )
% title('Displacement error at t_{end}')
ylabel('error')
xlabel('x')

E_displacement=sqrt(sum(volume_p(:,n_time_steps).*(u_X_p(:,n_time_steps)-...
    disp_X_exact(:,n_time_steps)).^2));
fprintf('L2-error in the u_x = %0.5e \n', E_displacement);

figure
hold on
plot( particles_X(1:n_particles_X,1),...
    v_X_p(1:n_particles_X,n_time_steps)-...
    velocity_X_func_exact(particles_X(1:n_particles_X,1),...
    particles_Y(1:n_particles_X,1),total_time) )
% title('x-velocity error')
ylabel('v')
xlabel('x')

E_vel=sqrt(sum(volume_p(:,n_time_steps).*(v_X_p(:,n_time_steps)-...
    velocity_X_func_exact(particles_X(:,1),...
    particles_Y(:,1),total_time)).^2));
fprintf('L2-error in the v_x = %0.5e \n', E_vel);

toc
profile viewer


