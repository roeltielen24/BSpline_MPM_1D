function [v_n] = initial_velocity(Xi,deg,n_dof,n_time_steps,velocity_func) 

% Determine initial velocity at degrees of freedom
if deg == 2  
    pos_p_glob2 = repmat([1/4 - (1/4)*(1/sqrt(3)),1/4 + (1/4)*(1/sqrt(3)),3/4 - (1/4)*(1/sqrt(3)), 3/4 + (1/4)*(1/sqrt(3))],1,nnz(Xi(2:end,1) - Xi(1:end-1)))'.*kron(nonzeros(Xi(2:end,1) - Xi(1:end-1)), ones(4,1))  + kron(Xi(1+deg:end-deg-1), ones(4,1));
    weight2 = 1/4*ones(length(pos_p_glob2),1);
    [N_vec2 B_vec2] = value_Bspline(Xi,deg,pos_p_glob2);
    weight3 = repmat([5/36, 8/36, 5/36, 5/36, 8/36, 5/36],1,nnz(Xi(2:end,1) - Xi(1:end-1)));
    pos_p_glob3 = repmat([  1/4 - (1/4)*(sqrt(3)/sqrt(5)),1/4,1/4 + (1/4)*(sqrt(3)/sqrt(5)), 3/4 - (1/4)*(sqrt(3)/sqrt(5)),3/4, 3/4 + (1/4)*(sqrt(3)/sqrt(5))],1,nnz(Xi(2:end,1) - Xi(1:end-1)))'.*kron(nonzeros(Xi(2:end,1) - Xi(1:end-1)), ones(6,1))  + kron(Xi(1+deg:end-deg-1), ones(6,1));   
    [N_vec3 B_vec3] = value_Bspline(Xi,deg,pos_p_glob3);
    for i = 1:n_dof
        f(i,1) = sum(weight2.*N_vec2(:,i).*velocity_func(pos_p_glob2));
        for j = 1:n_dof
            M1(i,j) = sum(weight3'.*N_vec3(:,i).*N_vec3(:,j));
        end
    end
end
if deg == 1   
        weight2 = 1/4*ones(4*(n_dof-deg),1);
        pos_p_glob2 = repmat([1/4 - (1/4)*(1/sqrt(3)),1/4 + (1/4)*(1/sqrt(3)),3/4 - (1/4)*(1/sqrt(3)), 3/4 + (1/4)*(1/sqrt(3))],1,nnz(Xi(2:end,1) - Xi(1:end-1)))'.*kron(nonzeros(Xi(2:end,1) - Xi(1:end-1)), ones(4,1))  + kron(Xi(1+deg:end-deg-1), ones(4,1));
        [N_vec2 B_vec2] = value_Bspline(Xi,deg,pos_p_glob2);
        weight1 = 1*ones(n_dof-deg,1);
        pos_p_glob1 =  repmat(1/2,1,n_dof-deg)'.*kron(nonzeros(Xi(2:end,1) - Xi(1:end-1)), ones(1,1))  + kron(Xi(1+deg:end-deg-1), ones(1,1)) ; 
        [N_vec, B_vec] = value_Bspline(Xi,deg,pos_p_glob1);
        for i = 1:n_dof
            f(i,1) = sum(weight1.*N_vec(:,i).*velocity_func(pos_p_glob1));
            M1(i,i) = sum(weight1.*N_vec(:,i));
        end
end

% Determine initial conditions
v_n = zeros(n_dof,n_time_steps-1);
v_n(:,1) = M1\f;
v_n(1,1) = 0;
v_n(end,1) = 0;
end