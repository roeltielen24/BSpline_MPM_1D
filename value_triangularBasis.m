%% Function value_triangularBasis
% Input: particle positions, the triangulation and the number of triangles
% Output: Value of material points in basis functions and value of the
% derivative to X and Y of the basis functions in the material points

function [N_vec, B_vec_X, B_vec_Y, nodes_active, triangles_active] = ...
    value_triangularBasis(particles_X,particles_Y,triangles)
%%
n_vertices=size(triangles.Points,1);
n_particles=numel(particles_X);

% N_vec  =zeros(n_particles,n_vertices);
% B_vec_X=zeros(n_particles,n_vertices);
% B_vec_Y=zeros(n_particles,n_vertices);

%List with particles with nonzero entries 
p_temp=repmat(1:n_particles,3,1);
p=reshape(p_temp,3*n_particles,1);

[p_tr,BC] = pointLocation(triangles,particles_X,particles_Y);
n=triangles.ConnectivityList(p_tr,:);
n=reshape(n',3*n_particles,1);
% p_tr=reshape(repmat(p_tr,1,3)',3*n_particles,1);

nodes_active = unique(n);
triangles_active = unique(p_tr);
%Check if any particles are outside the domain
if (sum(isnan(p_tr))~=0)
    disp('Particle outside domain')
end

N_vec = reshape(BC',3*n_particles,1);
N_vec = sparse(p,n,N_vec,n_particles,n_vertices);

V_X=triangles.Points(:,1);
V_Y=triangles.Points(:,2);
V1_X=V_X(triangles.ConnectivityList(p_tr,1));
V2_X=V_X(triangles.ConnectivityList(p_tr,2));
V3_X=V_X(triangles.ConnectivityList(p_tr,3));
V1_Y=V_Y(triangles.ConnectivityList(p_tr,1));
V2_Y=V_Y(triangles.ConnectivityList(p_tr,2));
V3_Y=V_Y(triangles.ConnectivityList(p_tr,3));

denom=( ( V2_Y-V3_Y ).*( V1_X-V3_X ) + ( V3_X-V2_X ).*( V1_Y-V3_Y ) );

deta1_dx=(V2_Y-V3_Y)./denom;
deta2_dx=(V3_Y-V1_Y)./denom;
deta3_dx=(V1_Y-V2_Y)./denom;

deta1_dy=(V3_X-V2_X)./denom;
deta2_dy=(V1_X-V3_X)./denom;
deta3_dy=(V2_X-V1_X)./denom;

B_vec_X=reshape([deta1_dx,deta2_dx,deta3_dx]',3*n_particles,1);
B_vec_Y=reshape([deta1_dy,deta2_dy,deta3_dy]',3*n_particles,1);

B_vec_X=sparse(p,n,B_vec_X,n_particles,n_vertices);
B_vec_Y=sparse(p,n,B_vec_Y,n_particles,n_vertices);

end
