%% Function value_triangularBasis
% Input: particle positions, the triangulation and the number of triangles
% Output: Value of material points in basis functions and value of the
% derivative to X and Y of the basis functions in the material points

function [N_vec, B_vec_X, B_vec_Y] = ...
    value_triangularBasis(particles_X,particles_Y,n_particles,triangles)
%%
n_vertices=size(triangles.Points,1);

N_vec  =zeros(n_particles,n_vertices);
B_vec_X=zeros(n_particles,n_vertices);
B_vec_Y=zeros(n_particles,n_vertices);

[point_loc,bary_coords] = pointLocation(triangles,particles_X,particles_Y);

for i=1:n_particles
    N_vec(i,triangles.ConnectivityList(point_loc(i),:))=bary_coords(i,:);
    
    %The barycentric coordinates expressed in x and y (source:
    %https://en.wikipedia.org/wiki/Barycentric_coordinate_system)
    %eta1=((y2-y3)(x-x3)+(x3-x2)(y-y3))/((y2-y3)(x1-x3)+(x3-x2)(y1-y3))
    %eta2=((y3-y1)(x-x3)+(x1-x3)(y-y3))/((y2-y3)(x1-x3)+(x3-x2)(y1-y3))
    %eta3=1-eta1-eta2
    %From these the derivatives to x and y can easily be derived
    
    %x and y locations of the vertices of the current triangle
    x=triangles.Points(triangles.ConnectivityList(point_loc(i),:),1);
    y=triangles.Points(triangles.ConnectivityList(point_loc(i),:),2);
    det=(y(2)-y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3));
    
    B_vec_X(i,triangles.ConnectivityList(point_loc(i),1))=...
        (y(2)-y(3))/det;
    B_vec_Y(i,triangles.ConnectivityList(point_loc(i),1))=...
        (x(3)-x(2))/det;
    
    B_vec_X(i,triangles.ConnectivityList(point_loc(i),2))=...
        (y(3)-y(1))/det;
    B_vec_Y(i,triangles.ConnectivityList(point_loc(i),2))=...
        (x(1)-x(3))/det;
    
    B_vec_X(i,triangles.ConnectivityList(point_loc(i),3))=...
        (y(1)-y(2))/det;
    B_vec_Y(i,triangles.ConnectivityList(point_loc(i),3))=...
        (x(2)-x(1))/det;
end
