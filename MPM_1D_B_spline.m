%% MPM_1D_B_spline Version 2-11-2016
%  Roel Tielen (Based on code of Lisa Wobbes), TU Delft
%  This file performs the MPM computation

function[u_n, v_p, v_n, mass_p,u_p,E_kin,E_pot,E_grav,E_trac,...
    stress_p,strain_p] = MPM_1D_B_spline(constant,flag,pos_p_glob,pos_p_loc, n_e, h, n_ep,...
     t_step, n_time_steps, total_time, E_kin,E_pot,E_grav,E_trac,...
     weight,Xi,deg,displacement_func, velocity_func,stress_func,n_p,Total_F) 
 
%% Mesh and knot element connections
n_dof = n_e + deg;                                                          % Number of degrees of freedom  
mesh  = Xi(1+deg:end-deg);                                                  % Knot span mesh
elements_nodes = [mesh(1:end-1) mesh(2:end)];                               % Store knot span position                                   
elements_nodes_index = [(1:n_e)'  (2:n_e+1)' ];                             % Store knot span index

elements_particles = zeros(n_e,n_p);                                        % Element where particle is located
for i=1:n_e
    for j = (i-1)*n_ep+1:i*n_ep
        elements_particles(i,j) = 1;                                        % Initial configuration
    end
end
 
% Initial particle properties
volume_p = (1/n_ep)*h*ones(n_p,1);                                          % Volume of particle
density_p = constant.density*ones(n_p,1);                                   % Density of particle
mass_p = volume_p.*density_p;                                               % Mass of particle
f_grav_p = mass_p*constant.g;                                               % Gravitational force particle
f_trac_p = zeros(n_p,1);                                                    % Traction force particle
f_trac_p(end) = constant.load;                                              % Application load

% Vectors to store information
stress_p = zeros(n_p,n_time_steps);                                         % Particle stress
strain_p = zeros(n_p,n_time_steps);                                         % Particle strain

v_n = initial_velocity(Xi,deg,n_dof,n_time_steps,velocity_func);            % Initial velocity at DOF
v_p = zeros(n_p, n_time_steps-1);                                           % Particle velocity                
v_p(:,1) = velocity_func(pos_p_glob);                                       % Initial particle velocity    

u_n = zeros(n_dof, n_time_steps-1);                                         % Displacement DOF
u_n(:,1) = zeros(n_dof,1);                                                  % Initial displacement at DOF
u_p = zeros(n_p, n_time_steps-1);                                           % Displacement particles
u_p(:,1) = zeros(n_p,1);                                                    % Initial displacement particles


        

%% Time integration
for s = 1:n_time_steps-1
    
    [N_vec, B_vec] = value_Bspline(Xi,deg,pos_p_glob);                      % Function values at particle positions
    
    if flag.splines == 1
        loc_gp = sort(repmat([1:2:3]/4,1,2)) ...                                    % Local position Gauss points
            + repmat([-1/sqrt(3) 1/sqrt(3)]/4,1,2);                         % ...
        pos_p_glob_gp = repmat(loc_gp,1,nnz(Xi(2:end,1) - Xi(1:end-1)))'... % Global position Gauss points
         .*kron(nonzeros(Xi(2:end,1) - Xi(1:end-1)), ones(length(loc_gp)... % ...
            ,1))  + kron(Xi(1+deg:end-deg-1), ones(length(loc_gp),1));      % ...
        if flag.both_ends_fixed == 0
            Xi_low  = max(Xi(Xi < pos_p_glob(end) + volume_p(end)/2));      % Knot before right boundary continuum
            Xi_high = min(Xi(Xi > Xi_low));                                 % Knot after right boundary continuum 
            if Xi_low == constant.height
                Xi_high = constant.height;
                pos_p_glob_gp = pos_p_glob_gp(pos_p_glob_gp < Xi_low);
            else
                pos_p_glob_gp = pos_p_glob_gp(pos_p_glob_gp < Xi_high);
                pos_p_glob_gp(end-3:end) = loc_gp*((pos_p_glob(end) + volume_p(end)/2)-Xi_low) + Xi_low;
            end
            
            % Use the correct integration weights
            Xi2 = Xi(Xi <= Xi_high);
            omega = repmat(nonzeros(Xi2(2:end) - Xi2(1:end-1)),4,1);
            omega(end-3:end) = pos_p_glob(end) + volume_p(end)/2 - Xi_low;
            omega = omega/4;
        else
            omega = repmat(h,length(pos_p_glob_gp),1)/4;
        end
        stress_gp  = spline(pos_p_glob,stress_p(:,s),pos_p_glob_gp);        
        [N_vec_gp B_vec_gp] = value_Bspline(Xi,deg,pos_p_glob_gp);
    end
    
    for i = 1:n_dof
        if flag.splines == 1
            M(i,i) = sum(weight.*N_vec(:,i).*density_p);
            F_int(i,1) = sum(omega.*B_vec_gp(:,i).*stress_gp);
            P(i,1) = sum(weight.*N_vec(:,i).*v_p(:,s).*density_p);  
            F_grav(i,1) = sum(mass_p.*N_vec(:,i).*constant.g);
           else
            if flag.ULFEM == 1 
                weight = mesh(2:end)-mesh(1:end-1);
                M(i,i) = sum(mass_p.*N_vec(:,i));
                P(i,1) = sum(mass_p.*N_vec(:,i).*v_p(:,s));
                F_grav(i,1) = sum(mass_p.*N_vec(:,i).*constant.g);
                F_int(i,1) = sum(weight.*B_vec(:,i).*stress_p(:,s));
            else
            M(i,i) = sum(weight.*N_vec(:,i).*density_p);
            P(i,1) = sum(weight.*N_vec(:,i).*v_p(:,s).*density_p);
            F_grav(i,1) = sum(weight.*N_vec(:,i).*density_p*constant.g);
            F_int(i,1) = sum(weight.*B_vec(:,i).*stress_p(:,s));
            end

        end
    end
    clear i
    % Determine total mass at DOF
    sum(sum(M))
    
    % Determine element vector for traction force
    F_trac = zeros(n_dof,1);  
    el = particle_element(elements_particles(:,end)); 
    if  flag.dynamic == 1
        F_trac(el,1)   = N_vec(end,el)*f_trac_p(end,s+1);
        F_trac(el+1,1) = N_vec(end,el+1)*f_trac_p(end,s+1);
    else 
        F_trac = zeros(n_dof,1);
        F_trac(end) = constant.load;
    end
    clear el
      
    % Determine damping force
    F_damp = - sign(v_n(:,s))*constant.alpha.*abs(F_trac + F_grav - F_int);
    
    % Determine total force
    F = F_trac + F_grav - F_int + F_damp;
    
    % Apply boundary condition on momentum
    if flag.both_ends_fixed == 1 
        P(1) = 0;
        P(n_dof) = 0;
    else 
        P(1) = 0;
    end
    
    % Find inactive nodes
    non_active_nodes = find(diag(M)==0);
   
    % Determine nodal velocity
    if flag.lagranian == 1
        if flag.momentum == 1
        for n = 1:n_dof
            if ~any(non_active_nodes == n)
                v_n(n,s) = P(n)/M(n,n);
            end
        end
        end
    end
    clear n 
    
    % Calculate the nodal acceleration
    a_n = zeros(n_dof,1);
    if flag.lumped == 1 
       for n = 1:n_dof      
            if ~any(non_active_nodes == n)
                a_n(n) = F(n)/M(n,n) ;
            end
        end
        clear n
    else
        a_n(:,1) = M\F;
    end
    
    % Apply boundsry condition on nodal acceleration
    if flag.both_ends_fixed == 1
        a_n(1,1) = 0;
        a_n(n_dof,1) = 0;
    else 
        a_n(1,1) = 0;
    end

    % Update particle velocity
    v_p(:,s+1) = v_p(:,s) + t_step*N_vec*a_n;
     
    % Calculate updated momentum in active elements
    Mv = N_vec'*(v_p(:,s+1).*mass_p);
    
    % Apply boundary condition on updated momentum
    if flag.both_ends_fixed == 1
        Mv(1)=0;
        Mv(n_dof) = 0;
    else
        Mv(1) = 0;
    end
    
    % Calculate the nodal velocities at new time step
    if flag.splines == 1
        velocity_gp2  = spline(pos_p_glob,v_p(:,s+1),pos_p_glob_gp);
        density_gp2  = spline(pos_p_glob,density_p,pos_p_glob_gp);
    end
      
    if flag.lagranian == 0
        for n = 1:n_dof
            if ~any(non_active_nodes == n)
               if flag.splines == 0
                    v_n(n,s+1) = Mv(n)/M(n,n);
               else
                    vtemp_n(n,s+1) = sum(omega.*N_vec_gp(:,n).*velocity_gp2.*density_gp2);
                    for j = 1:n_dof
                        M1(n,j) = sum(omega.*N_vec_gp(:,n).*N_vec_gp(:,j).*density_gp2);
                    end
                    
                end
            end
        end
        clear n
   
    if flag.splines == 1
        %Find DOF which are not zero
        nzm  = find(sum(M,2) ~= 0);
         
        % Now apply Richardson iteration!
        Pf = vtemp_n(nzm,s+1);
        v_n(nzm,s+1) = Pf(nzm)./sum(M1(nzm,nzm))';
        while (norm(Pf - M1(nzm,nzm)*v_n(nzm,s+1))/norm(v_n(nzm,s+1)) >= 1e-9)
            v_n(nzm,s+1) = (Pf + (diag(sum(M1(nzm,nzm))) - M1(nzm,nzm))*v_n(nzm,s+1))./sum(M1(nzm,nzm))';
        end
    end
        
    else
        v_n(:,s+1) = v_n(:,s) + a_n*t_step;
    end
    
    % Apply boundary condition on nodal velocity
    if flag.both_ends_fixed == 1
        v_n(1,s+1) = 0;
        v_n(n_dof,s+1) = 0;
    else
        v_n(1,s+1) = 0;
    end
     
    % Calculate the nodal incremental displacement
    du_n = v_n(:,s+1)*t_step;

    % Update the nodal displacement
    u_n(:,s+1) = u_n(:,s) + du_n;
    
    % Apply boundary condition on nodal displacement
    if flag.both_ends_fixed == 1
        u_n(1,s+1) = 0;
        u_n(n_dof,s+1) = 0;
    else
        u_n(1,s+1) = 0;
    end
    
    % Update knot vector according to obtained velocity - only deg = 1!
    if flag.ULFEM == 1
        [N_vec_up B_vec_up] = value_Bspline(Xi,deg,mesh);
        N_vec_up(end,end) = 1;
        mesh = mesh + N_vec_up*v_n(:,s+1)*t_step;
        Xi = [repmat(mesh(1),deg,1) ; mesh ; repmat(mesh(end),deg,1)];
       [N_vec, B_vec] = test_new_bspline(Xi,deg,pos_p_glob);
    end
    
    % Calculate the incremental particle strain
    dstrain_p = B_vec*du_n;
    strain_p(:,s+1)  = strain_p(:,s) + dstrain_p;
    
    % Update particle stress
    if flag.deformation == 1
        stress_p(:,s+1) = stress_p(:,s) + constant.E*dstrain_p - stress_p(:,s).*dstrain_p;
    else
        stress_p(:,s+1) = stress_p(:,s) + constant.E*dstrain_p;
    end
    
    % Update particle position and displacement
    u_p(:,s+1) = u_p(:,s) + N_vec*du_n;
    if flag.change_glob_pos == 1
        pos_p_glob = pos_p_glob + N_vec*du_n;
    end
 
   % Update the volume and density of the particles
   if flag.volume_update == 1
       volume_p  = (1 + dstrain_p).*volume_p;
       density_p = density_p./(1 + dstrain_p);
   end 
   
   % Update weights for integration (= volume)
   weight = volume_p; 
   
   % Check if there are particles that have left their original element
   if flag.change_glob_pos == 1
       for p = 1:n_p
           el = particle_element(elements_particles(:,p));
           if pos_p_glob(p) < elements_nodes(el,1)
               fprintf('Particle crossed lower border!')
               indicator(el,s) = 700;
               old_elem = el;
               if old_elem > 1
                   new_elem = old_elem - 1;
                   elements_particles(old_elem,p) = 0;
                   elements_particles(new_elem,p) = 1;
                   time_lower_crossing = s*t_step
               else
                   fprintf('Particle crossed the lower boundary!')
               end
               clear old_elem new_elem
           elseif pos_p_glob(p) > elements_nodes(el,2)
               fprintf('Particle crossed upper border!')
               indicator(el,s) = 700;
               old_elem = el;
               if old_elem < n_e
                   new_elem = old_elem + 1;
                   elements_particles(old_elem, p) = 0;
                   elements_particles(new_elem,p) = 1;
                   time_upper_crossing = s*t_step
               else
                   fprintf('Particle crossed the upper bundary!')
               end
               clear old_eleme new_elem
           end
       end
       clear p
   end
    
    % Change the local positions of the particles
    if flag.change_loc_pos == 1
        for p = 1:n_p
            pos_p_loc(p) = (pos_p_glob(p) - elements_nodes(particle_element...
            (elements_particles(:,p)),1))/h;
        end
    end
    clear p
    
    % Check whether an empty cells occurs and print time
    particles_in_element = sum(elements_particles,2);
    
    for el=1:n_e
       if particles_in_element(el) == 0
           disp('Empty cell!');
           time_empty_cell = s*t_step
           number_elements_in_particle = sum(elements_particles,2);
       end
    end
    clear el
    
    % Determine the kinetic and potential energy of the system
    for n = 1:n_dof
        if flag.lumped == 0 
            for m = 1:n_dof
                E_kin(1,s+1) = E_kin(1,s+1) + 0.5*M(m,n)*v_n(n,s+1)*v_n(n,s+1);
            end
        else
            E_kin(1,s+1) = E_kin(1,s+1) + 0.5*M(n,n)*v_n(n,s+1)*v_n(n,s+1);
        end
        E_pot(1,s+1) = E_pot(1,s+1) + -0.5*F_int(n)*u_n(n,s+1);
        E_trac(1,s+1) = E_trac(1,s+1) + F_trac(n)*u_n(n,s+1);
        E_grav(1,s+1) = E_grav(1,s+1) + F_grav(n)*u_n(n,s+1);
    end
end

