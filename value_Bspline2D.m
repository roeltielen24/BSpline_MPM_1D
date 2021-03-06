%% Function value_Bspline
% Input: Bezier ordinate locations, bezier ordinates, particle positions,
% the triangulation and the number of triangles
% Output: Value of material points in basis functions and value of the
% derivative to X and Y of the basis functions in the material points

function [N_vec, B_vec_X, B_vec_Y] = ...
    value_Bspline2D(B_o_X,B_o_Y,B_o,particles_X,particles_Y,triangles,...
    n_triangles,n_dof)

n_vertices=n_dof/3;

N_vec  =zeros(length(particles_X),3*n_vertices);
B_vec_X=zeros(length(particles_X),3*n_vertices);
B_vec_Y=zeros(length(particles_X),3*n_vertices);

%T_temp contains the connectivity list of the subtriangles, of which there
%are 6*n_triangles
T_temp=reshape(1:3*6*n_triangles,3,[])';

%B_o_n is reshaped such that it becomes a long vector. B_o_n_temp is sorted
%as follows: first we sort by main triangle, then by subtriangle and last
%by the vertices. So B_o_n(tr,subtr,Ver) goes to 
%(B_o_n_temp((tr-1)*6*3+(subtr-1)*3+Ver).
B_o_X_temp=zeros(3*6*n_triangles,1);
B_o_Y_temp=zeros(3*6*n_triangles,1);
for i=1:n_triangles
    for j=1:6
        B_o_X_temp((i-1)*3*6+(j-1)*3+(1:3))=...
            squeeze(B_o_X(i,j,[1,3,5]));
        B_o_Y_temp((i-1)*3*6+(j-1)*3+(1:3))=...
            squeeze(B_o_Y(i,j,[1,3,5]));
    end
end

TR=triangulation(T_temp,B_o_X_temp,B_o_Y_temp);

%The subtriangle in which each material point lies
point_loc = pointLocation(TR,[particles_X,particles_Y]);

%In case a point is exactly on the boundary between two triangles, it can
%happen that is assigned to neither of the triangles due to the used
%method. In this case, we move the particle a tad bit such that it will be
%assigned to a triangle.
point_loc_NaN=find(isnan(point_loc));
i=0;
while (~isempty(point_loc_NaN))
    i=i+1;
    point_loc(point_loc_NaN)=pointLocation(TR,...
        particles_X(point_loc_NaN)+(rand-1/2)*(1+10^i*eps), ...
        particles_Y(point_loc_NaN)+(rand-1/2)*(1+10^i*eps) );
    if (i==15)
        fprintf('Cannot find the triangle of particle %i ',...
            point_loc_NaN(1))
        return
    end
    point_loc_NaN=find(isnan(point_loc));
end
    
%We know in which subtriangle the point is, and there are 6 subtriangles
%per main triangle, therefore, we can find the main triangle by dividing by
%6 and rounding up, due to the numbering of the subtriangles.
point_tr = ceil(point_loc/6);
point_subtr = mod(point_loc-1,6)+1;

%Barycentric coordinates of the particles of the subtriangles they are in
p_bary_coords=...
    cartesianToBarycentric(TR,point_loc,[particles_X,particles_Y]);

%loop over all main vertices

%For B_o(_n), note that the corners 1,2 and 3 of the subtriangle are stored 
%at places 1,3 and 5 respectively!
for v=1:n_vertices
    %find for which particles the degrees of freedom of vertex v are
    %nonzero
    [tri_with_v,~]= find(triangles.ConnectivityList==v);
    p_near_v= find(ismember(point_tr,tri_with_v));
    
    %loop over all degrees of freedom at that vertex
    for q=1:3
        %specify the degree of freedom
        n=3*(v-1)+q;

        for p=p_near_v'
            N_vec(p,n)=...
                B_o(v,q,point_tr(p),point_subtr(p),1)*...
                    p_bary_coords(p,1)^2 +...
                B_o(v,q,point_tr(p),point_subtr(p),2)*...
                    2*p_bary_coords(p,1)*p_bary_coords(p,2) +...
                B_o(v,q,point_tr(p),point_subtr(p),3)*...
                    p_bary_coords(p,2)^2 +...
                B_o(v,q,point_tr(p),point_subtr(p),4)*...
                    2*p_bary_coords(p,2)*p_bary_coords(p,3) +...
                B_o(v,q,point_tr(p),point_subtr(p),5)*...
                    p_bary_coords(p,3)^2 +...
                B_o(v,q,point_tr(p),point_subtr(p),6)*...
                    2*p_bary_coords(p,3)*p_bary_coords(p,1) ;
            
            %The derivatives easily derived from the explicit formula for
            %eta(x,y), which is given on wikipedia (amongst others)
            %https://en.wikipedia.org/wiki/Barycentric_coordinate_system
            denom=( ( B_o_Y(point_tr(p),point_subtr(p),3)-...
                      B_o_Y(point_tr(p),point_subtr(p),5) )*...
                    ( B_o_X(point_tr(p),point_subtr(p),1)-...
                      B_o_X(point_tr(p),point_subtr(p),5) ) ...
                      + ...
                    ( B_o_X(point_tr(p),point_subtr(p),5)-...
                      B_o_X(point_tr(p),point_subtr(p),3) )*...
                    ( B_o_Y(point_tr(p),point_subtr(p),1)-...
                      B_o_Y(point_tr(p),point_subtr(p),5) ) );
            
            %For calculating the derivatives to x and y. detan_di denotes 
            %the partial derivative of eta_n to i (x or y) in the particle 
            %point  
            eta1    = p_bary_coords(p,1);
            deta1_dx=(B_o_Y(point_tr(p),point_subtr(p),3)-...
                      B_o_Y(point_tr(p),point_subtr(p),5))/denom;
            deta1_dy=(B_o_X(point_tr(p),point_subtr(p),5)-...
                      B_o_X(point_tr(p),point_subtr(p),3))/denom;
            eta2    = p_bary_coords(p,2);
            deta2_dx=(B_o_Y(point_tr(p),point_subtr(p),5)-...
                      B_o_Y(point_tr(p),point_subtr(p),1))/denom;
            deta2_dy=(B_o_X(point_tr(p),point_subtr(p),1)-...
                      B_o_X(point_tr(p),point_subtr(p),5))/denom;
            eta3    = p_bary_coords(p,3);
            deta3_dx=(B_o_Y(point_tr(p),point_subtr(p),1)-...
                      B_o_Y(point_tr(p),point_subtr(p),3))/denom;
            deta3_dy=(B_o_X(point_tr(p),point_subtr(p),3)-...
                      B_o_X(point_tr(p),point_subtr(p),1))/denom;
            
            %The derivative to x of the basis function in p
            B_vec_X(p,n)=...
                B_o(v,q,point_tr(p),point_subtr(p),1)*...
                    2*eta1*deta1_dx +...
                B_o(v,q,point_tr(p),point_subtr(p),2)*...
                    2*(eta2*deta1_dx + eta1*deta2_dx) +...   
                B_o(v,q,point_tr(p),point_subtr(p),3)*...
                    2*eta2*deta2_dx +...
                B_o(v,q,point_tr(p),point_subtr(p),4)*...
                    2*(eta3*deta2_dx + eta2*deta3_dx) +... 
                B_o(v,q,point_tr(p),point_subtr(p),5)*...
                    2*eta3*deta3_dx +...
                B_o(v,q,point_tr(p),point_subtr(p),6)*...
                    2*(eta1*deta3_dx + eta3*deta1_dx);
                
            %The derivative to y of the basis function in p
            B_vec_Y(p,n)=...
                B_o(v,q,point_tr(p),point_subtr(p),1)*...
                    2*eta1*deta1_dy +...
                B_o(v,q,point_tr(p),point_subtr(p),2)*...
                    2*(eta2*deta1_dy + eta1*deta2_dy) +...   
                B_o(v,q,point_tr(p),point_subtr(p),3)*...
                    2*eta2*deta2_dy +...
                B_o(v,q,point_tr(p),point_subtr(p),4)*...
                    2*(eta3*deta2_dy + eta2*deta3_dy) +... 
                B_o(v,q,point_tr(p),point_subtr(p),5)*...
                    2*eta3*deta3_dy +...
                B_o(v,q,point_tr(p),point_subtr(p),6)*...
                    2*(eta1*deta3_dy + eta3*deta1_dy);
                   
        end %p
        
    end %q
    
end %v
 
end
