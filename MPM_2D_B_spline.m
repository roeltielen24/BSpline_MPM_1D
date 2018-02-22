%% MPM_1D_B_spline Version 2-11-2016
%  Pascal de Koster(Based on code of Roel Tielen and Lisa Wobbes), TU Delft
%  This file performs the MPM computation

function[u_X_n, u_Y_n, particles_X, particles_Y, v_X_p, v_Y_p,...
    v_X_n, v_Y_n, mass_p,u_X_p, u_Y_p,...
    E_kin,E_pot,E_grav,E_trac,...
    stress_p,strain_p,N_vec, B_vec_X, B_vec_Y] = ...
    MPM_2D_B_spline(constant,flag,n_dof,n_triangles,...
    n_particles,particles_X,particles_Y,vertices_X,vertices_Y,triangles,...
    area_weight, velocity_X_func, velocity_Y_func,t_step,n_time_steps)

%% Mesh and initialisation

%If splines are used as basis functions, store them with the Bezier
%Ordinates, otherwise use standard Lagrange basis tent functions
if flag.spline_basis==1
    [PS_tri_X,PS_tri_Y,B_o_X,B_o_Y,B_o] = ...
        Basis_functions_Bezier_ordinates(vertices_X,vertices_Y,triangles);
end

elements_particles = pointLocation(triangles,...
    particles_X(:,1),particles_Y(:,1));                                     % Triangle where particle is located

% Initial particle properties
volume_p = area_weight;                                                     % Volume of particle
density_p = constant.density*ones(n_particles,1);                           % Density of particle
mass_p = volume_p.*density_p;                                               % Mass of particle
f_grav_p = mass_p*constant.g;                                               % Gravitational force particle
f_trac_p = zeros(n_particles,1);                                            % Traction force particle
f_trac_p(end) = constant.load;                                              % Application load

% Vectors to store information
%stress_p(i,j,p,t)=sigma_{i,j} of particle p at time T_t.
stress_p = zeros(constant.dim,constant.dim,n_particles,n_time_steps);       % Particle stress
strain_p = zeros(constant.dim,constant.dim,n_particles,n_time_steps);       % Particle strain

v_X_p = zeros(n_particles, n_time_steps-1);                                 % Particle velocity     
v_Y_p = zeros(n_particles, n_time_steps-1);                                 % Particle velocity 
v_X_p(:,1) = velocity_X_func(particles_X(:,1),particles_Y(:,1),0);          % Initial particle velocity    
v_Y_p(:,1) = velocity_Y_func(particles_X(:,1),particles_Y(:,1),0);          % Initial particle velocity    

v_X_n = zeros(n_dof, n_time_steps-1);
v_Y_n = zeros(n_dof, n_time_steps-1);

u_X_n = zeros(n_dof, n_time_steps-1);                                       % Displacement DOF
u_X_n(:,1) = zeros(n_dof,1);                                                % Initial displacement at DOF
u_Y_n = zeros(n_dof, n_time_steps-1);                                       % Displacement DOF
u_Y_n(:,1) = zeros(n_dof,1);                                                % Initial displacement at DOF

u_X_p = zeros(n_particles, n_time_steps-1);                                 % Displacement particles
u_X_p(:,1) = zeros(n_particles,1);                                          % Initial displacement particles
u_Y_p = zeros(n_particles, n_time_steps-1);                                 % Displacement particles
u_Y_p(:,1) = zeros(n_particles,1);                                          % Initial displacement particles

E_kin_X=zeros(1,n_time_steps);
E_kin_Y=zeros(1,n_time_steps);
E_kin  =zeros(1,n_time_steps);
E_pot_X=zeros(1,n_time_steps);
E_pot_Y=zeros(1,n_time_steps);
E_pot  =zeros(1,n_time_steps);
E_trac_X=zeros(1,n_time_steps);
E_trac_Y=zeros(1,n_time_steps);
E_trac  =zeros(1,n_time_steps);
E_grav=zeros(1,n_time_steps);

E_kin_X(1)=0.5*sum(mass_p.*v_X_p(:,1).^2);
E_kin_Y(1)=0.5*sum(mass_p.*v_Y_p(:,1).^2);

%degrees of freedom of the vertices on the left and right boundary
if flag.spline_basis==1
    %The degrees of freedom important for the left and right boundary are
    %those corresponding to the Powell-Sabin triangle points on those
    %boundaries (Isogeometric analysis with Powell–Sabin splines for
    %advection–diffusion–reaction problems, Proposition 1, p.134 and p.137)
    dof_left  =find(abs(reshape(PS_tri_X',[],1)                )<1e-14);
    dof_right =find(abs(reshape(PS_tri_X',[],1)-constant.length)<1e-14);
    dof_bottom=find(abs(reshape(PS_tri_Y',[],1)                )<1e-14);
    dof_top   =find(abs(reshape(PS_tri_Y',[],1)-constant.width )<1e-14);
else
    %The degrees of freedom correspond to the vertices on the left and
    %right boundary
    dof_left=find(vertices_X==0);
    dof_right=find(vertices_X==constant.length);
    dof_bottom=find(vertices_Y==0);
    dof_top=find(vertices_Y==constant.width);
end

%degrees of freedom at both bottom or top
dof_bt=union(dof_bottom,dof_top);
    
%Bulk modulus
K=constant.E/( 3*(1-2*constant.Poisson_ratio) );
%Shear modulus
G=constant.E/( 2*(1+  constant.Poisson_ratio) );

%% Time integration

for s=1:n_time_steps-1
    %N_vec(p,n) contains the value of basis function n in point p. Likewise
    %for the derivative of the basis function to X and Y in B_vec_X/Y.
    if flag.spline_basis==1
        [N_vec, B_vec_X, B_vec_Y]=value_Bspline2D(B_o_X,B_o_Y,B_o,...
            particles_X(:,s),particles_Y(:,s),triangles,n_triangles,n_dof);                      % Function values at particle positions
    else
        [N_vec, B_vec_X, B_vec_Y] = ...
            value_triangularBasis(particles_X(:,s),particles_Y(:,s),...
            n_particles,triangles);
    end
    
    %Calculate the Mass matrix as \int(phi_i*rho*phi_j)
    %This can be approximated as Sum_n [w_n*(phi_i*rho*phi_j)|_{x_n}]
    
    M=zeros(n_dof);
    
    %Calculate the mass matrix by 
    M=N_vec'*(((area_weight.*density_p)*ones(1,n_dof)).*N_vec);
    
%     %Here an alternative way for making the mass matrix, that is by first
%     %constructing a spline through the particle points and calculate the
%     %density mass of the basis functions in gauss points and then integrate
%     
%     
%     
%     
    
    %Force vectors in the degrees of freedom, see page 11 of Roel's thesis
    F_grav_X_n = N_vec'*(area_weight.*density_p*constant.g);
    F_grav_Y_n = zeros(n_dof,1);
    F_int_X_n = B_vec_X'*(area_weight.*squeeze(stress_p(1,1,:,s)))+...
              B_vec_Y'*(area_weight.*squeeze(stress_p(2,1,:,s)));
    F_int_Y_n = B_vec_X'*(area_weight.*squeeze(stress_p(1,2,:,s)))+...
              B_vec_Y'*(area_weight.*squeeze(stress_p(2,2,:,s)));
    
    %Traction force on particles. When present, this is only applied to
    %particles at the top of the domain.
    F_trac_X_n=zeros(n_dof,1);
    F_trac_Y_n=zeros(n_dof,1);
    %!!! Load is not applied yet !!!
    %The problem lies in the that we do not know which particles are at the
    %end of the domain and therefore do not know on which particles to
    %apply the load.
    
    % Determine damping force
    F_damp_X_n = - sign(v_X_n(:,s))*constant.alpha.*...
        abs(F_trac_X_n + F_grav_X_n - F_int_X_n);
    F_damp_Y_n = - sign(v_Y_n(:,s))*constant.alpha.*...
        abs(F_trac_Y_n + F_grav_Y_n - F_int_Y_n);
    
    F_X_n = F_trac_X_n + F_grav_X_n - F_int_X_n + F_damp_X_n;
    F_Y_n = F_trac_Y_n + F_grav_Y_n - F_int_Y_n + F_damp_Y_n;
    
    % Calculate the nodal acceleration with boundary conditions
    % Apply boundary condition on nodal acceleration  
    M_a=M;
    F_X_a=F_X_n;
    F_Y_a=F_Y_n;
    
    M_a(dof_left,:)=0;
    M_a(:,dof_left)=0;
    M_a(dof_left,dof_left)=eye(length(dof_left));
    F_X_a(dof_left)=0;
    F_Y_a(dof_left)=0;
    if flag.both_ends_fixed == 1
        M_a(dof_right,:)=0;
        M_a(:,dof_right)=0;
        M_a(dof_right,dof_right)=eye(length(dof_right));
        F_X_a(dof_right)=0;
        F_Y_a(dof_right)=0;  
    end
    
    %Top and bottom have a free slip BC for the x-direction
    M_ax=M_a;

    %Add top and bottom Dirichlet boundary condtion for y-direction
    M_ay=M_a;
    M_ay(dof_bt,:)=0;
    M_ay(:,dof_bt)=0;
    M_ay(dof_bt,dof_bt)=eye(length(dof_bt));
    F_Y_a(dof_bt)=0;
    
    a_X_n = M_ax\F_X_a;
    a_Y_n = M_ay\F_Y_a;

    % Update particle velocity
    v_X_p(:,s+1) = v_X_p(:,s) + t_step*(N_vec*a_X_n);
    v_Y_p(:,s+1) = v_Y_p(:,s) + t_step*(N_vec*a_Y_n);
    
    % Calculate updated momentum (mass*velocity P=M*v) in active elements
    P_X_n = N_vec'*(volume_p.*density_p.*v_X_p(:,s+1));
    P_Y_n = N_vec'*(volume_p.*density_p.*v_Y_p(:,s+1));

    %Project the updated velocity onto the active elements
    M_v=M;
    P_X_n_v=P_X_n;
    P_Y_n_v=P_Y_n;
    
    M_v(dof_left,:)=0;
    M_v(:,dof_left)=0;
    M_v(dof_left,dof_left)=eye(length(dof_left));
    P_X_n_v(dof_left)=0;
    P_Y_n_v(dof_left)=0;
    if flag.both_ends_fixed == 1
        M_v(dof_right,:)=0;
        M_v(:,dof_right)=0;
        M_v(dof_right,dof_right)=eye(length(dof_right));
        P_X_n_v(dof_right)=0;
        P_Y_n_v(dof_right)=0;  
    end
    
    %Top and bottom have a free slip BC for the x-direction
    M_vx=M_v;

    %Add top and bottom Dirichlet boundary condtion for y-direction
    M_vy=M_v;
    M_vy(dof_bt,:)=0;
    M_vy(:,dof_bt)=0;
    M_vy(dof_bt,dof_bt)=eye(length(dof_bt));
    P_Y_n_v(dof_bt)=0;
    
    v_X_n(:,s+1) = M_vx\P_X_n_v;
    v_Y_n(:,s+1) = M_vy\P_Y_n_v;
    

    
%     
%     %plot the projection of the velocity
%     [X,Y]=meshgrid(0:constant.length/100:constant.length,...
%         0:constant.length/20:constant.width);
%     [Z,~,~]=value_Bspline2D(B_o_X,B_o_Y,B_o,...
%             reshape(X,[],1),reshape(Y,[],1),triangles,n_triangles,n_dof);
%     Z=Z*
%         
%     scatter(X,Y,Z);
%     
%     
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For checking how the projection works out                            %
    if s==1%n_time_steps-1                                                %
    resol=20;                                                             %
    figure                                                                %
    [xq,yq] = meshgrid(0:constant.length/resol:constant.length,...        %
        0:constant.width/resol:constant.width);                           %
    F = scatteredInterpolant(particles_X(:,s),...                         %
        particles_Y(:,s),v_X_p(:,s),'nearest');                           %
    z1=F(xq,yq);                                                          %
    plot3(particles_X(:,s),particles_Y(:,s),v_X_p(:,s),'mo')              %
    hold on                                                               %
    mesh(xq,yq,z1)                                                        %
    title('Nearest Neighbor')                                             %
    legend('Sample Points','Interpolated Surface','Location','NorthWest') %
                                                                          %
    figure                                                                %
    if flag.spline_basis==1                                               %
        [N_vec2, B_vec2, ~]=value_Bspline2D(B_o_X,B_o_Y,B_o,...           %
            reshape(xq,[],1),reshape(yq,[],1),triangles,...               %
            n_triangles,n_dof);                                           %
    else                                                                  %
        [N_vec2, B_vec2, ~]=value_triangularBasis(...                     %
            reshape(xq,[],1),reshape(yq,[],1),numel(xq),triangles);       %
    end                                                                   %
    z2=reshape(N_vec2*v_X_n(:,s+1),size(xq));                             %
%     z2=reshape(N_vec2*[zeros(12,1);1;zeros(n_dof-1-12,1)],size(xq));    %
    mesh(xq,yq,z2)                                                        %
    vel_X_exact2=velocity_X_func(xq,yq,s*t_step);                         %
    L2_error=1/resol^2*sum(sum(abs(z2-vel_X_exact2)));                    %
    fprintf('The L2 error is %.1e ', L2_error);                           %
    end                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Displacement from initial position of the particles
    du_X_n = v_X_n(:,s+1)*t_step;
    du_Y_n = v_Y_n(:,s+1)*t_step;
    
    %Update the nodal displacement
    u_X_n(:,s+1) = u_X_n(:,s) + du_X_n;
    u_Y_n(:,s+1) = u_Y_n(:,s) + du_X_n;
    
    %Calculate new strain, Eq 2.2 in Roel's thesis
    dstrain_11_p = t_step*(B_vec_X*v_X_n(:,s+1));
    dstrain_12_p = t_step*1/2*(B_vec_Y*v_X_n(:,s+1)+B_vec_X*v_Y_n(:,s+1));
    dstrain_21_p = dstrain_12_p;
    dstrain_22_p = t_step*(B_vec_Y*v_Y_n(:,s+1));
    
    strain_p(1,1,:,s+1)  = squeeze(strain_p(1,1,:,s)) + dstrain_11_p;
    strain_p(1,2,:,s+1)  = squeeze(strain_p(1,2,:,s)) + dstrain_12_p;
    strain_p(2,1,:,s+1)  = squeeze(strain_p(2,1,:,s)) + dstrain_21_p;
    strain_p(2,2,:,s+1)  = squeeze(strain_p(2,2,:,s)) + dstrain_22_p;
    
    %Calculate new stress. Eq 2.5 in Roel's thesis and 3.12 of MPM for
    %Geomechanical Problems
    %For 2D
    
    %First define the change due to Kirchhoff/Hill stress
    dstress_11_H_p=(2*G+(K-2/3*G))*dstrain_11_p+...
                   (    (K-2/3*G))*dstrain_22_p;
    dstress_22_H_p=(    (K-2/3*G))*dstrain_11_p+...
                   (2*G+(K-2/3*G))*dstrain_22_p;
    dstress_12_H_p=(2*G          )*dstrain_12_p;
    dstress_21_H_p=dstress_12_H_p;             
    
    %Calculate spin tensor, the antisymmetric part of the velocity gradient
    %tensor (eq 3.13) of MPM for Geomechanical Problems
    w_spin_11=zeros(n_particles,1);
    w_spin_12=t_step*1/2*(B_vec_Y*v_X_n(:,s+1)-B_vec_X*v_Y_n(:,s+1));
    w_spin_21=-w_spin_12;
    w_spin_22=zeros(n_particles,1);

    dstress_11_p=dstress_11_H_p ...
        -(dstrain_11_p+dstrain_22_p).*squeeze(stress_p(1,1,:,s)) ...
        + ( w_spin_11.*squeeze(stress_p(1,1,:,s)) + ...
            w_spin_12.*squeeze(stress_p(2,1,:,s)) ) ...
        - ( squeeze(stress_p(1,1,:,s)).*w_spin_11 + ...
            squeeze(stress_p(1,2,:,s)).*w_spin_21 );
    dstress_12_p=dstress_12_H_p ...
        -(dstrain_11_p+dstrain_22_p).*squeeze(stress_p(1,2,:,s)) ...
        + ( w_spin_11.*squeeze(stress_p(1,2,:,s)) + ...
            w_spin_12.*squeeze(stress_p(2,2,:,s)) ) ...
        - ( squeeze(stress_p(1,1,:,s)).*w_spin_12 + ...
            squeeze(stress_p(1,2,:,s)).*w_spin_22 );
    dstress_21_p=dstress_21_H_p ...
        -(dstrain_11_p+dstrain_22_p).*squeeze(stress_p(2,1,:,s)) ...
        + ( w_spin_21.*squeeze(stress_p(1,1,:,s)) + ...
            w_spin_22.*squeeze(stress_p(2,1,:,s)) ) ...
        - ( squeeze(stress_p(2,1,:,s)).*w_spin_11 + ...
            squeeze(stress_p(2,2,:,s)).*w_spin_21 );
    dstress_22_p=dstress_22_H_p ...
        -(dstrain_11_p+dstrain_22_p).*squeeze(stress_p(2,2,:,s)) ...
        + ( w_spin_21.*squeeze(stress_p(1,2,:,s)) + ...
            w_spin_22.*squeeze(stress_p(2,2,:,s)) ) ...
        - ( squeeze(stress_p(2,1,:,s)).*w_spin_12 + ...
            squeeze(stress_p(2,2,:,s)).*w_spin_22 );
    
    stress_p(1,1,:,s+1) = squeeze(stress_p(1,1,:,s)) + dstress_11_p;
    stress_p(1,2,:,s+1) = squeeze(stress_p(1,2,:,s)) + dstress_12_p;
    stress_p(2,1,:,s+1) = squeeze(stress_p(2,1,:,s)) + dstress_21_p;
    stress_p(2,2,:,s+1) = squeeze(stress_p(2,2,:,s)) + dstress_22_p;
    
    
    %VRAAG ROEL OVER PAGINA 14, TWEEDE VERGELIJKING (OVER STRESS)
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %For 3D, strain to stress, see wikipedia, Hooke's law for a
%     %conversion matrix
%     dstress_11_H_p=(2*G+(K-2/3*G))*dstrain_11_p+ ...
%                  (    (K-2/3*G))*dstrain_22_p+ ...
%                  (    (K-2/3*G))*dstrain_33_p;
%     dstress_22_H_p=(    (K-2/3*G))*dstrain_11_p+ ...
%                  (2*G+(K-2/3*G))*dstrain_22_p+ ...
%                  (    (K-2/3*G))*dstrain_33_p;
%     dstress_33_H_p=(    (K-2/3*G))*dstrain_11_p+ ...
%                  (    (K-2/3*G))*dstrain_22_p+ ...
%                  (2*G+(K-2/3*G))*dstrain_33_p;
%     dstress_12_H_p=(2*G          )*dstrain_12_p;
%     dstress_13_H_p=(2*G          )*dstrain_13_p;
%     dstress_23_H_p=(2*G          )*dstrain_23_p;
%     dstress_21_H_p=dstress_12_H_p;   
%     dstress_31_H_p=dstress_13_H_p;   
%     dstress_32_H_p=dstress_23_H_p
%     
%     %Calculate spin tensor, the antisymmetric part of the velocity gradient
%     %tensor (eq 3.13) of MPM for Geomechanical Problems
%     w_spin_11=zeros(n_particles,1);
%     w_spin_12=t_step*1/2*(B_vec_Y*v_X_n(:,s+1)-B_vec_X*v_Y_n(:,s+1));
%     w_spin_13=t_step*1/2*(B_vec_Z*v_X_n(:,s+1)-B_vec_X*v_Z_n(:,s+1));
%     w_spin_21=-w_spin_12;
%     w_spin_22=zeros(n_particles,1);
%     w_spin_23=t_step*1/2*(B_vec_Z*v_Y_n(:,s+1)-B_vec_Y*v_Z_n(:,s+1));
%     w_spin_31=-w_spin_13;
%     w_spin_32=-w_spin_23;
%     w_spin_33=zeros(n_particles,1);
%     
%     dstress_11_p=dstress_11_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(1,1,:,s)) ...
%         + ( w_spin_11.*squeeze(stress_p(1,1,:,s)) + ...
%             w_spin_12.*squeeze(stress_p(2,1,:,s)) + ...
%             w_spin_13.*squeeze(stress_p(3,1,:,s)) ) ...
%         - ( squeeze(stress_p(1,1,:,s)).*w_spin_11 + ...
%             squeeze(stress_p(1,2,:,s)).*w_spin_21 + ...
%             squeeze(stress_p(1,3,:,s)).*w_spin_31 );
%     dstress_12_p=dstress_12_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(1,2,:,s)) ...
%         + ( w_spin_11.*squeeze(stress_p(1,2,:,s)) + ...
%             w_spin_12.*squeeze(stress_p(2,2,:,s)) + ...
%             w_spin_13.*squeeze(stress_p(3,2,:,s)) ) ...
%         - ( squeeze(stress_p(1,1,:,s)).*w_spin_12 + ...
%             squeeze(stress_p(1,2,:,s)).*w_spin_22 + ...
%             squeeze(stress_p(1,3,:,s)).*w_spin_32 );
%     dstress_13_p=dstress_13_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(1,3,:,s)) ...
%         + ( w_spin_11.*squeeze(stress_p(1,3,:,s)) + ...
%             w_spin_12.*squeeze(stress_p(2,3,:,s)) + ...
%             w_spin_13.*squeeze(stress_p(3,3,:,s)) ) ...
%         - ( squeeze(stress_p(1,1,:,s)).*w_spin_13 + ...
%             squeeze(stress_p(1,2,:,s)).*w_spin_23 + ...
%             squeeze(stress_p(1,3,:,s)).*w_spin_33 );
%     dstress_21_p=dstress_21_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(2,1,:,s)) ...
%         + ( w_spin_21.*squeeze(stress_p(1,1,:,s)) + ...
%             w_spin_22.*squeeze(stress_p(2,1,:,s)) + ...
%             w_spin_23.*squeeze(stress_p(3,1,:,s)) ) ...
%         - ( squeeze(stress_p(2,1,:,s)).*w_spin_11 + ...
%             squeeze(stress_p(2,2,:,s)).*w_spin_21 + ...
%             squeeze(stress_p(2,3,:,s)).*w_spin_31 );
%     dstress_22_p=dstress_22_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(2,2,:,s)) ...
%         + ( w_spin_21.*squeeze(stress_p(1,2,:,s)) + ...
%             w_spin_22.*squeeze(stress_p(2,2,:,s)) + ...
%             w_spin_23.*squeeze(stress_p(3,2,:,s)) ) ...
%         - ( squeeze(stress_p(2,1,:,s)).*w_spin_12 + ...
%             squeeze(stress_p(2,2,:,s)).*w_spin_22 + ...
%             squeeze(stress_p(2,3,:,s)).*w_spin_32 );
%     dstress_23_p=dstress_23_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(2,3,:,s)) ...
%         + ( w_spin_21.*squeeze(stress_p(1,3,:,s)) + ...
%             w_spin_22.*squeeze(stress_p(2,3,:,s)) + ...
%             w_spin_23.*squeeze(stress_p(3,3,:,s)) ) ...
%         - ( squeeze(stress_p(2,1,:,s)).*w_spin_13 + ...
%             squeeze(stress_p(2,2,:,s)).*w_spin_23 + ...
%             squeeze(stress_p(2,3,:,s)).*w_spin_33 );        
%     dstress_31_p=dstress_31_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(3,1,:,s)) ...
%         + ( w_spin_31.*squeeze(stress_p(1,1,:,s)) + ...
%             w_spin_32.*squeeze(stress_p(2,1,:,s)) + ...
%             w_spin_33.*squeeze(stress_p(3,1,:,s)) ) ...
%         - ( squeeze(stress_p(3,1,:,s)).*w_spin_11 + ...
%             squeeze(stress_p(3,2,:,s)).*w_spin_21 + ...
%             squeeze(stress_p(3,3,:,s)).*w_spin_31 );
%     dstress_32_p=dstress_32_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(3,2,:,s)) ...
%         + ( w_spin_31.*squeeze(stress_p(1,2,:,s)) + ...
%             w_spin_32.*squeeze(stress_p(2,2,:,s)) + ...
%             w_spin_33.*squeeze(stress_p(3,2,:,s)) ) ...
%         - ( squeeze(stress_p(3,1,:,s)).*w_spin_12 + ...
%             squeeze(stress_p(3,2,:,s)).*w_spin_22 + ...
%             squeeze(stress_p(3,3,:,s)).*w_spin_32 );
%     dstress_33_p=dstress_33_H_p ...
%         -(dstrain_11_p+dstrain_22_p+dstrain_33_p)...
%             .*squeeze(stress_p(3,3,:,s)) ...
%         + ( w_spin_31.*squeeze(stress_p(1,3,:,s)) + ...
%             w_spin_32.*squeeze(stress_p(2,3,:,s)) + ...
%             w_spin_33.*squeeze(stress_p(3,3,:,s)) ) ...
%         - ( squeeze(stress_p(3,1,:,s)).*w_spin_13 + ...
%             squeeze(stress_p(3,2,:,s)).*w_spin_23 + ...
%             squeeze(stress_p(3,3,:,s)).*w_spin_33 );
% 
%     stress_p(1,1,:,s+1) = stress_p(1,1,:,s) + dstress_11_p;
%     stress_p(1,2,:,s+1) = stress_p(1,2,:,s) + dstress_12_p;
%     stress_p(1,3,:,s+1) = stress_p(1,3,:,s) + dstress_13_p;
%     stress_p(2,1,:,s+1) = stress_p(2,1,:,s) + dstress_21_p;
%     stress_p(2,2,:,s+1) = stress_p(2,2,:,s) + dstress_22_p;    
%     stress_p(2,3,:,s+1) = stress_p(2,3,:,s) + dstress_23_p;    
%     stress_p(3,1,:,s+1) = stress_p(3,1,:,s) + dstress_31_p;
%     stress_p(3,2,:,s+1) = stress_p(3,2,:,s) + dstress_32_p;
%     stress_p(3,3,:,s+1) = stress_p(3,3,:,s) + dstress_33_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update particle displacement
    u_X_p(:,s+1) = u_X_p(:,s) + N_vec*du_X_n;
    u_Y_p(:,s+1) = u_Y_p(:,s) + N_vec*du_Y_n;
    
    %Update particle positions
    particles_X(:,s+1) = particles_X(:,s) + N_vec*du_X_n;
    particles_Y(:,s+1) = particles_Y(:,s) + N_vec*du_Y_n;
    
    % Update the volume and density of the particles
    volume_p  = (1 + dstrain_11_p+dstrain_22_p).*volume_p;
    density_p = density_p./(1 + dstrain_11_p + dstrain_22_p);
    
    % Update weights for integration (= volume)
    area_weight = volume_p;
    
    % Determine the kinetic and potential energy of the system
    E_kin_X(1,s+1)=0.5*(sum(M,1)*v_X_n(:,s+1).^2);
    E_kin_Y(1,s+1)=0.5*(sum(M,1)*v_Y_n(:,s+1).^2);
    E_kin(1,s+1)=E_kin_X(1,s+1)+E_kin_Y(1,s+1);
    
    E_pot_X(1,s+1)=-0.5*F_int_X_n'*u_X_n(:,s+1);
    E_pot_Y(1,s+1)=-0.5*F_int_Y_n'*u_Y_n(:,s+1);
    E_pot(1,s+1)=E_pot_X(1,s+1)+E_pot_Y(1,s+1);
    
    E_trac_X(1,s+1) = F_trac_X_n'*u_X_n(:,s+1);
    E_trac_Y(1,s+1) = F_trac_Y_n'*u_Y_n(:,s+1);
    E_trac(1,s+1) = E_trac_X(1,s+1)+E_trac_Y(1,s+1);
    
    E_grav(1,s+1) = E_grav(1,s+1) + F_grav_X_n'*u_X_n(:,s+1);    
end



















end
