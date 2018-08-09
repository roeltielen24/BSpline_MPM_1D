%% Function value_Bspline
% Input: Bezier ordinate locations, bezier ordinates, particle positions,
% the triangulation and the number of triangles Output: Value of material
% points in basis functions and value of the derivative to X and Y of the
% basis functions in the material points

function [N_vec, B_vec_X, B_vec_Y, nodes_active,triangles_active] = ...
    value_Bspline2D(B_o_X,B_o_Y,B_o,...
    particles_X,particles_Y,triangles,n_triangles,n_dof)

n_particles=length(particles_X);

%List with particles with nonzero entries 
p_temp=repmat(1:n_particles,9,1);
p=reshape(p_temp,9*n_particles,1);

%Find the nonzero basis functions in the particle locations. For each
%particle, there are 9 basis function that are non-zero in particle p
%(each particle lies in a triangle with 3 vertices, of whic each has 3
%basis functions)
p_tr = pointLocation(triangles,[particles_X,particles_Y]);
p_st = zeros(n_particles,1);
BC = zeros(n_particles,3);

%Check if any particles are outside the domain
if (sum(isnan(p_tr))~=0)
    disp('Particle outside domain')
end

n=triangles.ConnectivityList(p_tr,:);
n_v=[n(:,1),n(:,1),n(:,1),n(:,2),n(:,2),n(:,2),n(:,3),n(:,3),n(:,3)];
n_v=reshape(n_v',9*n_particles,1);
n_q=repmat([1;2;3],3*n_particles,1);
n=[3*n(:,1)-2,3*n(:,1)-1,3*n(:,1),...
   3*n(:,2)-2,3*n(:,2)-1,3*n(:,2),...
   3*n(:,3)-2,3*n(:,3)-1,3*n(:,3)];
n=reshape(n',9*n_particles,1);

nodes_active = unique(n);
triangles_active = unique(p_tr);
%Next find the subtriangle each particle is in
for tr = 1:n_triangles
    %The triangle is subdivided into subtriangles, numbered as follows
    %       6  
    %     / | \
    %    7. | .5
    %   /_.'1'._\
    %  2----3----4
    %Here, node number 2, 4 and 6 are respectively the nodes 1, 2 and 3 in
    %the main triangle
    CL=[1 1 1 1 1 1;...
        2 3 4 5 6 7;...
        3 4 5 6 7 2]';
    Points = [B_o_X(tr,1,1) , B_o_Y(tr,1,1) ;...
              B_o_X(tr,:,3)', B_o_Y(tr,:,3)' ]; 
    SubTri = triangulation(CL,Points);
    
    %Find the subtriangles and the respective barycentric coordinates of
    %each particle in triangle tr
    p_tr_i=find(p_tr == tr); 
    if isempty(p_tr_i) %then triangle tr does not contain any points
        continue
    end 
    [p_st_temp,BC_temp] = pointLocation(SubTri,...
        particles_X(p_tr_i),particles_Y(p_tr_i));
    
    NaNi = find(isnan(p_st_temp));
    if ~isempty(NaNi)
        for i=1:10 
            %In case a particle's triangle cannot be determined, it must be
            %on an edge. Pull it a bit towards the center vertex to make
            %sure it is inside a triangle.
            [sti2,BC2] =pointLocation(SubTri,...
                particles_X(p_tr_i(NaNi)) + ...
                100^i*eps*(Points(1,1)-particles_X(p_tr_i(NaNi))) ,...
                particles_Y(p_tr_i(NaNi)) + ...
                100^i*eps*(Points(1,2)-particles_Y(p_tr_i(NaNi))) );
                p_st_temp(NaNi)=sti2;
            BC_temp(NaNi,:)=BC2;
            NaNi = find(isnan(p_st_temp));
            if isempty(NaNi)
                break
            end
        end
    end   
    p_st(p_tr_i) = p_st_temp;
    BC(p_tr_i,:)=BC_temp;
end

if ~isempty(NaNi)
    disp(['Subtriangle index not found, probably there is a particle ',...
        'on an internal subtriangle edge'])
end

BC=repmat(BC,1,9);
BC=reshape(BC',3,9*n_particles)';

p_tr=reshape(repmat(p_tr,1,9)',9*n_particles,1);
p_st=reshape(repmat(p_st,1,9)',9*n_particles,1);

ind1=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,  ones(9*n_particles,1));
ind2=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,2*ones(9*n_particles,1));
ind3=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,3*ones(9*n_particles,1));
ind4=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,4*ones(9*n_particles,1));
ind5=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,5*ones(9*n_particles,1));
ind6=sub2ind(size(B_o),n_v,n_q,p_tr,p_st,6*ones(9*n_particles,1));

N_vec =...
    B_o(ind1).*BC(:,1).^2 +...
    B_o(ind2).*2.*BC(:,1).*BC(:,2) +...
    B_o(ind3).*BC(:,2).^2 +...
    B_o(ind4).*2.*BC(:,2).*BC(:,3) +...
    B_o(ind5).*BC(:,3).^2 +...
    B_o(ind6).*2.*BC(:,3).*BC(:,1);

N_vec = sparse(p,n,N_vec,n_particles,n_dof);

ind_der1=sub2ind(size(B_o_Y),p_tr,p_st,  ones(9*n_particles,1));
ind_der3=sub2ind(size(B_o_Y),p_tr,p_st,3*ones(9*n_particles,1));
ind_der5=sub2ind(size(B_o_Y),p_tr,p_st,5*ones(9*n_particles,1));

%The derivatives easily derived from the explicit formula for
%eta(x,y), which is given on wikipedia (amongst others)
%https://en.wikipedia.org/wiki/Barycentric_coordinate_system
denom=( ( B_o_Y(ind_der3)-B_o_Y(ind_der5) ).*...
        ( B_o_X(ind_der1)-B_o_X(ind_der5) ) ...
          + ...
        ( B_o_X(ind_der5)-B_o_X(ind_der3) ).*...
        ( B_o_Y(ind_der1)-B_o_Y(ind_der5) ) );

%For calculating the derivatives to x and y. detan_di denotes 
%the partial derivative of eta_n to i (x or y) in the particle 
%point  
eta1    = BC(:,1);
deta1_dx=(B_o_Y(ind_der3)-B_o_Y(ind_der5))./denom;
deta1_dy=(B_o_X(ind_der5)-B_o_X(ind_der3))./denom;
eta2    = BC(:,2);
deta2_dx=(B_o_Y(ind_der5)-B_o_Y(ind_der1))./denom;
deta2_dy=(B_o_X(ind_der1)-B_o_X(ind_der5))./denom;
eta3    = BC(:,3);
deta3_dx=(B_o_Y(ind_der1)-B_o_Y(ind_der3))./denom;
deta3_dy=(B_o_X(ind_der3)-B_o_X(ind_der1))./denom;
    

%The derivative to x of the basis function in p
B_vec_X=...
    B_o(ind1).*2.*eta1.*deta1_dx +...
    B_o(ind2).*2.*(eta2.*deta1_dx + eta1.*deta2_dx) +...   
    B_o(ind3).*2.*eta2.*deta2_dx +...
    B_o(ind4).*2.*(eta3.*deta2_dx + eta2.*deta3_dx) +... 
    B_o(ind5).*2.*eta3.*deta3_dx +...
    B_o(ind6).*2.*(eta1.*deta3_dx + eta3.*deta1_dx);

%The derivative to y of the basis function in p
B_vec_Y=...
    B_o(ind1).*2.* eta1.*deta1_dy +...
    B_o(ind2).*2.*(eta2.*deta1_dy + eta1.*deta2_dy) +...   
    B_o(ind3).*2.* eta2.*deta2_dy +...
    B_o(ind4).*2.*(eta3.*deta2_dy + eta2.*deta3_dy) +... 
    B_o(ind5).*2.* eta3.*deta3_dy +...
    B_o(ind6).*2.*(eta1.*deta3_dy + eta3.*deta1_dy);

B_vec_X = sparse(p,n,B_vec_X,n_particles,n_dof);
B_vec_Y = sparse(p,n,B_vec_Y,n_particles,n_dof);

end




