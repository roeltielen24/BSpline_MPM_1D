%% MPM_1D_B_spline Version 2-11-2016
%  Pascal de Koster(Based on code of Roel Tielen and Lisa Wobbes), TU Delft
%  This file performs the MPM computation

function[u_X_n, u_Y_n, particles_X, particles_Y, v_X_p, v_Y_p,...
    v_X_n, v_Y_n, mass_p,u_X_p, u_Y_p,...
    E_kin,E_pot,E_grav,E_trac,...
    stress11_p,stress12_p,stress21_p,stress22_p, ...
    strain11_p,strain12_p,strain21_p,strain22_p, ...
    N_vec, B_vec_X, B_vec_Y, volume_p] = ...
    MPM_2D_B_spline(constant,flag,n_dof,n_triangles,...
    n_particles,particles_X,particles_Y,vertices_X,vertices_Y,triangles,...
    area_weight, velocity_X_func, velocity_Y_func,t_step,n_time_steps)

%% Mesh and initialisation

%If splines are used as basis functions, store them with the Bezier
%Ordinates, otherwise use standard Lagrange basis tent functions
if flag.spline_basis==1
    [PS_tri_X,PS_tri_Y,B_o_X,B_o_Y,B_o] = ...
        Basis_functions_Bezier_ordinates(vertices_X,vertices_Y,...
        triangles,flag);
end

% figure(1)
% scatter(Q_X(v,:),Q_Y(v,:),'filled')

%Find Gauss points in the refined grid in case of Gauss Points
if flag.Gauss_integration ==1
    if flag.spline_basis ==1
        %Find the Gauss points in the refined triangulation. The refined
        %triangulation still has to be constructed though, so this is done
        %first.
        %Each subtriangle has the vertices with Bezier coordinate numbers
        %1,3 and 5  (2, 4 and 6 are halfway on the edges)
        ST_PointsX=permute(B_o_X(:,:,[1,3,5]),[3,2,1]);
        ST_PointsX=reshape(ST_PointsX,[],1);
        ST_PointsY=permute(B_o_Y(:,:,[1,3,5]),[3,2,1]);
        ST_PointsY=reshape(ST_PointsY,[],1);
        ST_ConnectivityList=( [1:3:3*6*n_triangles;...
                               2:3:3*6*n_triangles;...
                               3:3:3*6*n_triangles ]' );
        SubTriangles=triangulation(ST_ConnectivityList,...
            ST_PointsX,ST_PointsY);
        [GaussPoints_X,GaussPoints_Y,Gauss_weight]=...
            Gauss_Points(SubTriangles);
        GaussPoints_X=reshape(GaussPoints_X',[],1);
        GaussPoints_Y=reshape(GaussPoints_Y',[],1);
        Gauss_weight=reshape(Gauss_weight',[],1);
    else
        %Find the Gauss points in the original triangulation
        [GaussPoints_X,GaussPoints_Y,Gauss_weight]=...
        Gauss_Points(triangles);
        GaussPoints_X=reshape(GaussPoints_X',[],1);
        GaussPoints_Y=reshape(GaussPoints_Y',[],1);
        Gauss_weight=reshape(Gauss_weight',[],1);
    end
end


if flag.Gauss_integration == 1
    %Plot Gauss points
    hold on
    scatter(GaussPoints_Y,GaussPoints_X)
end

% Initial particle properties
volume_p=zeros(n_particles,n_time_steps);
volume_p(:,1) = area_weight;                                                % Volume of particle
density_p = constant.density*ones(n_particles,1);                           % Density of particle
mass_p = volume_p(:,1).*density_p;                                          % Mass of particle

% Vectors to store information
%stress_p(i,j,p,t)=sigma_{i,j} of particle p at time T_t.
% Particle stress
stress11_p = zeros(n_particles,n_time_steps);
stress12_p = zeros(n_particles,n_time_steps);
stress21_p = zeros(n_particles,n_time_steps);
stress22_p = zeros(n_particles,n_time_steps);
% Particle strain
strain11_p = zeros(n_particles,n_time_steps);
strain12_p = zeros(n_particles,n_time_steps);
strain21_p = zeros(n_particles,n_time_steps);
strain22_p = zeros(n_particles,n_time_steps);

v_X_p = zeros(n_particles, n_time_steps-1);                                 % Particle velocity     
v_Y_p = zeros(n_particles, n_time_steps-1);                                 % Particle velocity 
v_X_p(:,1) = velocity_X_func(particles_X(:,1),particles_Y(:,1));            % Initial particle velocity    
v_Y_p(:,1) = velocity_Y_func(particles_X(:,1),particles_Y(:,1));            % Initial particle velocity    

v_X_n = zeros(n_dof, n_time_steps);
v_Y_n = zeros(n_dof, n_time_steps);

a_X_n = zeros(n_dof, 1);
a_Y_n = zeros(n_dof, 1);

u_X_n = zeros(n_dof, n_time_steps);                                         % Displacement DOF
u_X_n(:,1) = zeros(n_dof,1);                                                % Initial displacement at DOF
u_Y_n = zeros(n_dof, n_time_steps);                                         % Displacement DOF
u_Y_n(:,1) = zeros(n_dof,1);                                                % Initial displacement at DOF

u_X_p = zeros(n_particles, n_time_steps);                                   % Displacement particles
u_X_p(:,1) = zeros(n_particles,1);                                          % Initial displacement particles
u_Y_p = zeros(n_particles, n_time_steps);                                   % Displacement particles
u_Y_p(:,1) = zeros(n_particles,1);                                          % Initial displacement particles

E_kin_X=zeros(1,n_time_steps);
E_kin_X(1,1)=sum(density_p.*area_weight.*v_X_p(:,1));
E_kin_Y=zeros(1,n_time_steps);
E_kin  =zeros(1,n_time_steps);
E_kin(1,1)=E_kin_X(1,1)+E_kin_Y(1,1);
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


%In case of Gauss integration, the integration points do not move, and the
%values of the basis functions only have to be calculated once
if flag.Gauss_integration==1
    %N_vec_gp(gp,n) contains the value of basis function n in Gauss point gp.
    %Likewise for the derivative of the basis function to X and Y in
    %B_vec_X/Y.
    if flag.spline_basis==1
        [N_vec_gp, B_vec_X_gp, B_vec_Y_gp, ~ ]=value_Bspline2D(...
            B_o_X,B_o_Y,B_o,GaussPoints_X,GaussPoints_Y,...
            triangles,n_triangles,n_dof); 
    else
        [N_vec_gp, B_vec_X_gp, B_vec_Y_gp, ~ ] = ...
            value_triangularBasis(GaussPoints_X,GaussPoints_Y,triangles);
    end
end

%% Time integration

for s=1:n_time_steps-1
    fprintf('s = %i \n', s)
%     if s == 269
%     end
%     if s == 174
%         s
%     end
%     if s==130
%     end
%     if s == 100
%     end
%     if s == 1377
%     end
%     if s == 9%890
%     end
    
    %N_vec(p,n) contains the value of basis function n in point p. Likewise
    %for the derivative of the basis function to X and Y in B_vec_X/Y.
    if flag.spline_basis==1
        [N_vec, B_vec_X, B_vec_Y, nodes_active,triangles_active]=...
            value_Bspline2D(...
            B_o_X,B_o_Y,B_o,...
            particles_X(:,s),particles_Y(:,s),triangles,n_triangles,...
            n_dof); 
    else
        [N_vec, B_vec_X, B_vec_Y, nodes_active,triangles_active] = ...
            value_triangularBasis(particles_X(:,s),particles_Y(:,s),...
            triangles);
    end
    
    if flag.DisableNodes == 1
        M_cons = N_vec'*(((area_weight.*density_p)*ones(1,n_dof)).*N_vec);
        nodes_zero_mass=find(diag(M_cons<=1e-2));
        %Filter out nodes that have hardly any mass
        nodes_active=setdiff(nodes_active,nodes_zero_mass);
    end
    
%     if ~any(nodes_active==55)
%         disp('node 55 not active')
%     end
%     if ~any(nodes_active==56)
%         disp('node 56 not active')
%     end
%     if ~any(nodes_active==57)
%         disp('node 57 not active')
%     end
    
    %Find the active vertices and triangles
    vertices_active = unique(ceil(nodes_active/3));
    vertices_inactive = setdiff(1:length(vertices_X),vertices_active);
    triangles_inactive = setdiff(1:n_triangles,triangles_active);
    vertices_active_boundary = setdiff(...
        unique(triangles.ConnectivityList(triangles_inactive,:)),...
        vertices_inactive);
    vertices_active_boundary = reshape(vertices_active_boundary,[],1);
    triangles_at_boundary = find(sum(ismember(...
                                         triangles.ConnectivityList,...
                                         vertices_active_boundary),2));
    triangles_active_boundary = setdiff(triangles_at_boundary,...
                                        triangles_inactive);
    particles_in_boundary_triangles = find(ismember(...
        pointLocation(triangles,particles_X(:,s),particles_Y(:,s)),...
        triangles_active_boundary));
    
    if s==114
    end
    %Find nodes from vertices
    if flag.spline_basis == 1   %PS-splines has 3 nodes per vertex
        if ~isempty(vertices_active_boundary)
            nodes_active_boundary=3*repmat(vertices_active_boundary',[3,1])...
                + [1;2;3]*ones(1,length(vertices_active_boundary)) -3;
            nodes_active_boundary = reshape(nodes_active_boundary,[],1);
        end
    else    %For Lagrange MPM, the vertices correspond to the nodes
        nodes_active_boundary = vertices_active_boundary;
    end
    
    
    %Calculate the Mass matrix as \int(phi_i*rho*phi_j)
    %This can be approximated as Sum_n [w_n*(phi_i*rho*phi_j)|_{x_n}]
    %Calculate the mass matrix by 
    if flag.lumped == 0     %Consistent system
        M = N_vec'*(((area_weight.*density_p)*ones(1,n_dof)).*N_vec);
    elseif flag.lumped == 2     %Partial lumping
        M = N_vec'*(((area_weight.*density_p)*ones(1,n_dof)).*N_vec);
        M(nodes_active_boundary,:)=0;
        M(nodes_active_boundary,nodes_active_boundary) = ...
            diag(N_vec(:,nodes_active_boundary)'*(area_weight.*density_p));
    else    %Lumping
        M = diag(N_vec'*(area_weight.*density_p));
        M_cons = N_vec'*(((area_weight.*density_p)*ones(1,n_dof)).*N_vec);
        
%         figure; surf(M); title('M')
    end        
    
    %Gravitational force in the degrees of freedom
    F_grav_X_n = N_vec'*(area_weight.*density_p*constant.g);
    F_grav_Y_n = zeros(n_dof,1);  
    
    %Stress tensors in the degrees of freedom, see page 11 of Roel's thesis
    if flag.Gauss_integration == 1
        if flag.TLS == 1    %Use Taylor Least Squares reconstruction
            conservation = 0;   %No conservation for sigma
            stress11_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                stress11_p(:,s),volume_p,GaussPoints_X,GaussPoints_Y,...
                Gauss_weight,triangles_active_boundary, ...
                triangles_inactive, triangles,conservation);
            stress12_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                stress12_p(:,s),volume_p,GaussPoints_X,GaussPoints_Y,...
                Gauss_weight,triangles_active_boundary, ...
                triangles_inactive, triangles,conservation);
            stress21_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                stress21_p(:,s),volume_p,GaussPoints_X,GaussPoints_Y,...
                Gauss_weight,triangles_active_boundary, ...
                triangles_inactive, triangles,conservation);
            stress22_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                stress22_p(:,s),volume_p,GaussPoints_X,GaussPoints_Y,...
                Gauss_weight,triangles_active_boundary, ...
                triangles_inactive, triangles,conservation);
            
            %In the inactive elements and the boundary elements, the values
            %in the gauss points are set to zero. In the boundary elements,
            %the integration is done with the original particles
            B_vec_X_temp = [B_vec_X_gp;...
                            B_vec_X(particles_in_boundary_triangles,:)];
            B_vec_Y_temp = [B_vec_Y_gp;...
                            B_vec_Y(particles_in_boundary_triangles,:)];
            weight_temp = [Gauss_weight;...
                           volume_p(particles_in_boundary_triangles)];
            stress11_temp = [stress11_gp;...
                             stress11_p(particles_in_boundary_triangles,s)];
            stress12_temp = [stress12_gp;...
                             stress12_p(particles_in_boundary_triangles,s)];
            stress21_temp = [stress21_gp;...
                             stress21_p(particles_in_boundary_triangles,s)];
            stress22_temp = [stress22_gp;...
                             stress22_p(particles_in_boundary_triangles,s)];

            F_int_X_n = B_vec_X_temp'*(weight_temp.*stress11_temp)+...
                        B_vec_Y_temp'*(weight_temp.*stress21_temp);
            F_int_Y_n = B_vec_X_temp'*(weight_temp.*stress12_temp)+...
                        B_vec_Y_temp'*(weight_temp.*stress22_temp);
        else    %Use cubic spline reconstruction
            f_stress11=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),stress11_p(:,s),'linear','nearest');
            stress11_gp = f_stress11(GaussPoints_X,GaussPoints_Y);
            f_stress12=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),stress12_p(:,s),'linear','nearest');
            stress12_gp = f_stress12(GaussPoints_X,GaussPoints_Y);
            f_stress21=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),stress21_p(:,s),'linear','nearest');
            stress21_gp = f_stress21(GaussPoints_X,GaussPoints_Y);
            f_stress22=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),stress22_p(:,s),'linear','nearest');
            stress22_gp = f_stress22(GaussPoints_X,GaussPoints_Y);
            
            F_int_X_n = B_vec_X_gp'*(Gauss_weight.*stress11_gp)+...
                        B_vec_Y_gp'*(Gauss_weight.*stress21_gp);
            F_int_Y_n = B_vec_X_gp'*(Gauss_weight.*stress12_gp)+...
                        B_vec_Y_gp'*(Gauss_weight.*stress22_gp);
        end
    else
        F_int_X_n = B_vec_X'*(area_weight.*stress11_p(:,s))+...
                    B_vec_Y'*(area_weight.*stress21_p(:,s));
        F_int_Y_n = B_vec_X'*(area_weight.*stress12_p(:,s))+...
                    B_vec_Y'*(area_weight.*stress22_p(:,s));
    end
    
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
%     M_a=M_cons;

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
    
    %Add top and bottom Dirichlet boundary condtion for y-acceleration
    M_ay=M_a;
    M_ay(dof_bt,:)=0;
    M_ay(:,dof_bt)=0;
    M_ay(dof_bt,dof_bt)=eye(length(dof_bt));
    F_Y_a(dof_bt)=0;
    
    %Remove inactive elements
    M_ax=M_ax(nodes_active,nodes_active);
    F_X_a=F_X_a(nodes_active);
    M_ay=M_ay(nodes_active,nodes_active);
    F_Y_a=F_Y_a(nodes_active);
    
%     if numel(nodes_active)~=3*length(vertices_X)
%         disp('There is an inactive node')
%         fprintf('at s = %i \n', s)
%     end
    
    %Solve for the acceleration in the degrees of freedom    
    if flag.lumped == 0     %Consistent system
        a_X_n(nodes_active) = M_ax\F_X_a;
        a_Y_n(nodes_active) = M_ay\F_Y_a;
    elseif flag.lumped == 2     %Boundary nodes are lumped
        a_X_n(nodes_active) = M_ax\F_X_a;
        a_Y_n(nodes_active) = M_ay\F_Y_a;
    else    %Normal lumping
        a_X_n(nodes_active) = F_X_a./diag(M_ax);
        a_Y_n(nodes_active) = F_Y_a./diag(M_ay);
        
%         a_X_n(nodes_active) = M_ax\F_X_a;
%         a_Y_n(nodes_active) = M_ay\F_Y_a;
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot the projection of the acceleration
    if mod(s,ceil(n_time_steps/10))==1
    [X,Y]=meshgrid((0.01:99.01)*constant.length/100,...
        (1:19)*constant.width/20);
    if (flag.spline_basis==1)
        [Z,~,~,~]=value_Bspline2D(B_o_X,B_o_Y,B_o,...
            reshape(X,[],1),reshape(Y,[],1),triangles,n_triangles,n_dof);
    else
        [Z, ~ , ~,~] = ...
            value_triangularBasis(reshape(X,[],1),reshape(Y,[],1),triangles);
    end
%     Z=reshape(Z*a_X_n,size(X));
    Z=reshape(Z*a_X_n,size(X));
    
    k=ceil(s/round(n_time_steps/10));
    figure(9)
    subplot(2,5,k)
    surf(X,Y,Z,'EdgeColor','none');
    title(sprintf('a_n, t = %.3f', s*t_step))
%     axis([0 constant.length 0 constant.width -1 1])
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update particle velocity
    v_X_p(:,s+1) = v_X_p(:,s) + t_step*(N_vec*a_X_n);
    v_Y_p(:,s+1) = v_Y_p(:,s) + t_step*(N_vec*a_Y_n);
    
    % Calculate updated momentum (mass*velocity P=M*v) in active elements
    if flag.Gauss_integration == 1
        if flag.TLS == 1
            conservation = 0;
            if flag.Cubic_TLS == 0
                %Quadratic
                P_X_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                    v_X_p(:,s+1).*density_p,volume_p,GaussPoints_X,...
                    GaussPoints_Y,Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
                P_Y_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                    v_Y_p(:,s+1).*density_p,volume_p,GaussPoints_X,...
                    GaussPoints_Y,Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
                density_gp=interpolate_TLS(particles_X(:,s),particles_Y(:,s),...
                    density_p,volume_p,GaussPoints_X,GaussPoints_Y,...
                    Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
            else
                %Cubic
                P_X_gp=interpolate_TLS_cubic(particles_X(:,s),particles_Y(:,s),...
                    v_X_p(:,s+1).*density_p,volume_p,GaussPoints_X,...
                    GaussPoints_Y,Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
                P_Y_gp=interpolate_TLS_cubic(particles_X(:,s),particles_Y(:,s),...
                    v_Y_p(:,s+1).*density_p,volume_p,GaussPoints_X,...
                    GaussPoints_Y,Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
                density_gp=interpolate_TLS_cubic(particles_X(:,s),particles_Y(:,s),...
                    density_p,volume_p,GaussPoints_X,GaussPoints_Y,...
                    Gauss_weight,triangles_active_boundary, ...
                    triangles_inactive, triangles,conservation);
            end
            
            %In the inactive elements and the boundary elements, the values
            %in the gauss points are set to zero. In the boundary elements,
            %the integration is done with the original particles
            N_vec_temp = [N_vec_gp;...
                          N_vec(particles_in_boundary_triangles,:)];
            weight_temp = [Gauss_weight;...
                           volume_p(particles_in_boundary_triangles)];
            P_X_temp = [P_X_gp;...
                        v_X_p(particles_in_boundary_triangles,s+1).*...
                        density_p(particles_in_boundary_triangles)];
            P_Y_temp = [P_Y_gp;...
                        v_Y_p(particles_in_boundary_triangles,s+1).*...
                        density_p(particles_in_boundary_triangles)];
            density_temp = [density_gp;...
                            density_p(particles_in_boundary_triangles)];
            F_int_X_n = B_vec_X_temp'*(weight_temp.*stress11_temp)+...
                        B_vec_Y_temp'*(weight_temp.*stress21_temp);
            F_int_Y_n = B_vec_X_temp'*(weight_temp.*stress12_temp)+...
                        B_vec_Y_temp'*(weight_temp.*stress22_temp);
            
            P_X_n = N_vec_temp'*(weight_temp.*P_X_temp);
            P_Y_n = N_vec_temp'*(weight_temp.*P_Y_temp);

            M_gp = N_vec_temp'*(((weight_temp.*density_temp)*...
                                ones(1,n_dof)).*N_vec_temp);  
            
        else 
%             f_v_X=scatteredInterpolant(particles_X(:,s),...
%                 particles_Y(:,s),v_X_p(:,s+1),'nearest','nearest');
%             v_X_gp = f_v_X(GaussPoints_X,GaussPoints_Y);
%             f_v_Y=scatteredInterpolant(particles_X(:,s),...
%                 particles_Y(:,s),v_Y_p(:,s+1),'nearest','nearest');
%             v_Y_gp = f_v_Y(GaussPoints_X,GaussPoints_Y);
            f_P_X=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),v_X_p(:,s+1).*density_p,...
                'linear','nearest');
            P_X_gp = f_P_X(GaussPoints_X,GaussPoints_Y);
            f_P_Y=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),v_Y_p(:,s+1).*density_p,...
                'linear','nearest');
            P_Y_gp = f_P_Y(GaussPoints_X,GaussPoints_Y);
            
            f_density=scatteredInterpolant(particles_X(:,s),...
                particles_Y(:,s),density_p,'linear','nearest');
            density_gp = f_density(GaussPoints_X,GaussPoints_Y);
            
            P_X_n = N_vec_gp'*(Gauss_weight.*P_X_gp);
            P_Y_n = N_vec_gp'*(Gauss_weight.*P_Y_gp);

            M_gp = N_vec_gp'*(((Gauss_weight.*density_gp)*...
                                ones(1,n_dof)).*N_vec_gp);  
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %plot the projection of the density
%         if mod(s,ceil(n_time_steps/10))==1
%         figure(12)
%         subplot(2,5,k)
%         surf(X,Y,Z,'EdgeColor','none');
%         title(sprintf('\rho_n, t = %.3f', s*t_step))
%         scatter3(GaussPoints_X,GaussPoints_Y,density_gp,100,...
%                 density_gp,'filled');
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        P_X_n = N_vec'*(volume_p(:,s).*density_p.*v_X_p(:,s+1));
        P_Y_n = N_vec'*(volume_p(:,s).*density_p.*v_Y_p(:,s+1));
    end
    
    %Project the updated velocity onto the active elements
    %First, make a copy of the mass matrix and momentum vectors for adding
    %boundary conditions
    if flag.Gauss_integration == 1
        M_v = M_gp;
    else
        M_v=M;
%         M_v = M_cons;
    end
    
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
    
    %Remove inactive elements
    M_vx=M_vx(nodes_active,nodes_active);
    P_X_n_v=P_X_n_v(nodes_active);
    M_vy=M_vy(nodes_active,nodes_active);
    P_Y_n_v=P_Y_n_v(nodes_active);
        
    %Solve for the new nodal velocities
    if flag.lumped == 0    %Consistent system
        v_X_n(nodes_active,s+1) = M_vx\P_X_n_v;
        v_Y_n(nodes_active,s+1) = M_vy\P_Y_n_v;
    elseif flag.lumped == 2     %Boundary nodes are lumped
        v_X_n(nodes_active,s+1) = M_vx\P_X_n_v;
        v_Y_n(nodes_active,s+1) = M_vy\P_Y_n_v;
    else     %Normal lumping
        v_X_n(nodes_active,s+1) = P_X_n_v./diag(M_vx);
        v_Y_n(nodes_active,s+1) = P_Y_n_v./diag(M_vy);
%         v_X_n(nodes_active,s+1) = M_vx\P_X_n_v;
%         v_Y_n(nodes_active,s+1) = M_vy\P_Y_n_v;
    end

%     %Compare consistent and lumped velocity in the nodes    
%     v_X_n_cons = M_vx\P_X_n_v;
    v_X_n_lump = P_X_n_v./diag(M_ax);
%     
%     if mod(s,ceil(n_time_steps/10))==1
%         figure;
%         hold on
%         plot(v_X_n_cons)
%         plot(v_X_n_lump)
%     end
    
%     %Compare consistent and lumped acceleration in the nodes
%     a_X_n_cons = M_vx\F_X_a;
%     a_X_n_lump = M_ax\F_X_a;
%     
%     figure;
%     hold on
%     plot(a_X_n_cons)
%     plot(a_X_n_lump)
%     hold off
%     legend('consistent','lumped')
%     disp(s)
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot the projection of the velocity
    if mod(s,ceil(n_time_steps/10))==1
    [X,Y]=meshgrid((0.01:99.01)*constant.length/100,...
        (1:19)*constant.width/20);
    if (flag.spline_basis==1)
        [Z,~,~,~]=value_Bspline2D(B_o_X,B_o_Y,B_o,...
            reshape(X,[],1),reshape(Y,[],1),triangles,n_triangles,n_dof);
    else
        [Z, ~ , ~, ~] = ...
            value_triangularBasis(reshape(X,[],1),reshape(Y,[],1),triangles);
    end
%     Z=reshape(Z*v_X_n(:,s+1),size(X));
    Z=reshape(Z(:,nodes_active)*v_X_n_lump,size(X));
    
    k=ceil(s/round(n_time_steps/10));
    figure(10)
    subplot(2,5,k)
    surf(X,Y,Z,'EdgeColor','none');
    title(sprintf('v_n, t = %.3f', s*t_step))
%     axis([0 constant.length 0 constant.width -1 1])
    disp('')
    
    figure(11)
    subplot(2,5,k)
    scatter3(particles_X(:,s),particles_Y(:,s),v_X_p(:,s),100,...
        v_X_p(:,s),'filled')
    title(sprintf('v_p, t = %.3f', s*t_step))
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Plot the velocity projection and calculate the L2-error              %
%     if s==1%n_time_steps-1                                                %
%     resol=20;                                                             %
%     figure                                                                %
%     [xq,yq] = meshgrid(0:constant.length/resol:constant.length,...        %
%         0:constant.width/resol:constant.width);                           %
%     F = scatteredInterpolant(particles_X(:,s),...                         %
%         particles_Y(:,s),v_X_p(:,s),'nearest');                           %
%     z1=F(xq,yq);                                                          %
%     plot3(particles_X(:,s),particles_Y(:,s),v_X_p(:,s),'mo')              %
%     hold on                                                               %
%     mesh(xq,yq,z1)                                                        %
%     title('Nearest Neighbor')                                             %
%     legend('Sample Points','Interpolated Surface','Location','NorthWest') %
%                                                                           %
%     figure                                                                %
%     if flag.spline_basis==1                                               %
%         [N_vec2, ~, ~,~]=value_Bspline2D(B_o_X,B_o_Y,B_o,...                %
%             reshape(xq,[],1),reshape(yq,[],1),triangles,...               %
%             n_triangles,n_dof);                                           %
%     else                                                                  %
%         [N_vec2, ~, ~, ~]=value_triangularBasis(...                          %
%             reshape(xq,[],1),reshape(yq,[],1),numel(xq),triangles);       %
%     end                                                                   %
%     z2=reshape(N_vec2*v_X_n(:,s+1),size(xq));                             %
% %     z2=reshape(N_vec2*[zeros(12,1);1;zeros(n_dof-1-12,1)],size(xq));    %
%     mesh(xq,yq,z2)                                                        %
%     vel_X_exact2=velocity_X_func(xq,yq,s*t_step);                         %
%     L2_error=1/resol^2*sum(sum(abs(z2-vel_X_exact2)));                    %
%     fprintf('The L2 error is %.1e ', L2_error);                           %
%     end                                                                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Displacement from initial position of the particles
    du_X_n = v_X_n(:,s+1)*t_step;
    du_Y_n = v_Y_n(:,s+1)*t_step;
    
    %Update the nodal displacement
    u_X_n(:,s+1) = u_X_n(:,s) + du_X_n;
    u_Y_n(:,s+1) = u_Y_n(:,s) + du_X_n;
    
    %Calculate new strain, Eq 2.2 in Roel's thesis
    dstrain11_p = t_step*(B_vec_X*v_X_n(:,s+1));
    dstrain12_p = t_step*1/2*(B_vec_Y*v_X_n(:,s+1)+...
                              B_vec_X*v_Y_n(:,s+1));
    dstrain21_p = dstrain12_p;
    dstrain22_p = t_step*(B_vec_Y*v_Y_n(:,s+1));
    
    strain11_p(:,s+1)  = strain11_p(:,s) + dstrain11_p;
    strain12_p(:,s+1)  = strain12_p(:,s) + dstrain12_p;
    strain21_p(:,s+1)  = strain21_p(:,s) + dstrain21_p;
    strain22_p(:,s+1)  = strain22_p(:,s) + dstrain22_p;
    
    %Calculate new stress. Eq 2.5 in Roel's thesis and 3.12 of MPM for
    %Geomechanical Problems
    %For 2D
    
    %First define the change due to Kirchhoff/Hill stress
    dstress11_H_p=(2*G+(K-2/3*G))*dstrain11_p+...
                  (    (K-2/3*G))*dstrain22_p;
    dstress22_H_p=(    (K-2/3*G))*dstrain11_p+...
                  (2*G+(K-2/3*G))*dstrain22_p;
    dstress12_H_p=(2*G          )*dstrain12_p;
    dstress21_H_p=dstress12_H_p;             
    
    %Calculate spin tensor, the antisymmetric part of the velocity gradient
    %tensor (eq 3.13) of MPM for Geomechanical Problems
    w_spin_11=sparse(n_particles,1);
    w_spin_12=t_step*1/2*(B_vec_Y*v_X_n(:,s+1)-...
                          B_vec_X*v_Y_n(:,s+1) );
    w_spin_21=-w_spin_12;
    w_spin_22=sparse(n_particles,1);

    dstress11_p=dstress11_H_p ...
        -(dstrain11_p+dstrain22_p).*stress11_p(:,s) ...
        + ( w_spin_11.*stress11_p(:,s) + ...
            w_spin_12.*stress21_p(:,s) ) ...
        - ( stress11_p(:,s).*w_spin_11 + ...
            stress12_p(:,s).*w_spin_21 );
    dstress12_p=dstress12_H_p ...
        -(dstrain11_p+dstrain22_p).*stress12_p(:,s) ...
        + ( w_spin_11.*stress12_p(:,s) + ...
            w_spin_12.*stress22_p(:,s) ) ...
        - ( stress11_p(:,s).*w_spin_12 + ...
            stress12_p(:,s).*w_spin_22 );
    dstress21_p=dstress21_H_p ...
        -(dstrain11_p+dstrain22_p).*stress21_p(:,s) ...
        + ( w_spin_21.*stress11_p(:,s) + ...
            w_spin_22.*stress21_p(:,s) ) ...
        - ( stress21_p(:,s).*w_spin_11 + ...
            stress22_p(:,s).*w_spin_21 );
    dstress22_p=dstress22_H_p ...
        -(dstrain11_p+dstrain22_p).*stress22_p(:,s) ...
        + ( w_spin_21.*stress12_p(:,s) + ...
            w_spin_22.*stress22_p(:,s) ) ...
        - ( stress21_p(:,s).*w_spin_12 + ...
            stress22_p(:,s).*w_spin_22 );
        
    stress11_p(:,s+1) = stress11_p(:,s) + dstress11_p;
    stress12_p(:,s+1) = stress12_p(:,s) + dstress12_p;
    stress21_p(:,s+1) = stress21_p(:,s) + dstress21_p;
    stress22_p(:,s+1) = stress22_p(:,s) + dstress22_p;
    
    % Update particle displacement
    du_X_p = N_vec*du_X_n;
    du_Y_p = N_vec*du_Y_n;
    
    u_X_p(:,s+1) = u_X_p(:,s) + du_X_p;
    u_Y_p(:,s+1) = u_Y_p(:,s) + du_Y_p;
        
    %Update particle positions
    particles_X(:,s+1) = particles_X(:,s) + du_X_p;
    particles_Y(:,s+1) = particles_Y(:,s) + du_Y_p;
    
    % Update the volume and density of the particles
    volume_p(:,s+1)  = (1 + dstrain11_p+dstrain22_p).*volume_p(:,s);
    density_p = density_p./(1 + dstrain11_p + dstrain22_p);
    
    % Update weights for integration (= volume)
    area_weight = volume_p(:,s+1);
    
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

% %Check the interpolation to the Gauss points
% if flag.Gauss_integration == 1
%     figure;
%     hold on;  
%     scatter3(particles_X(:,s),particles_Y(:,s),stress11_p(:,s));
%     scatter3(GaussPoints_X,GaussPoints_Y,stress11_gp,'r')
% end



end
