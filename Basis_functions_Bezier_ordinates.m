%% Function value_Bspline
% Input: X and Y locations vertices of tiangles (V_X,V_Y) and the triangles
% (Tri) Output: The X and Y locations of  the Powell-Sabin triangles (Q_X,
% Q_Y), the 3D matrix with the X and Y locations of the Bezier Ordinates
% (B_o_X, B_o_Y). These are defined as B_o_n(triangle,i,j), where n can be
% X or Y, and B_o_n denotes the location of the Bezier ordinate in main
% triangle 'triangle', of which subtriangle 'i' (counted counterclockwise,
% starting with the subtriangle with the following three vertices: the
% first vertex of the main triangle, the center vertex of the main triangle
% and the vertex on the edge between the first and second main triangle
% vertices. Then j denotes which location at in this subtriangle, counting
% counterclockwise, always starting at the location at the center vertex of
% the main triangle.
% Finally, B_o is a 5D matrix B_O(V,p,tri,i,j), denoting the
% Bezier ordinates of the basis function of main vertex V, of which basis
% function p (1,2 or 3), at the location of main triangle tri, and i and j
% defining subtriangle i, of which location j, identically to B_O_n.

function [Q_X,Q_Y,B_o_X,B_o_Y,B_o] = ...
    Basis_functions_Bezier_ordinates(V_X,V_Y,triangles,flag)

%% Construct Powell Sabin refinement

%Allow figures to be plotted
fig=0;

Num_Tri=size(triangles,1);
Tri_nb = neighbors(triangles);
Num_Ver=length(V_X);

%Generate points in the center of each triangle with incenter, which is the
%intersection point of the three bisectors of the angles of the triangle. 
%http://mathworld.wolfram.com/Incenter.html
if flag.grid_regularity==1
    V_c_X=mean(V_X(triangles.ConnectivityList),2);
    V_c_Y=mean(V_Y(triangles.ConnectivityList),2);
    V_c=[V_c_X,V_c_Y];
else
    V_c=incenter(triangles);
    V_c_X=V_c(:,1);
    V_c_Y=V_c(:,2);
end
if fig==1
    figure(1)
    hold on
    scatter(V_c(:,1),V_c(:,2))
    hold off
end

%Define V_edge_X/Y(Tri,Ver) as the vertex in Triangle Tri on the OPPOSITE
%side of the Vertex it belongs to. This way, the numbering is also
%unambiguous for tetrahedrons.
V_edge_X=zeros(Num_Tri,3);
V_edge_Y=zeros(Num_Tri,3);

%Determine the inner triangle points on the edges of the main triangle.
for i=1:Num_Tri
    for j=1:3   %Number of edges
        if isnan(Tri_nb(i,j)) %No neighbouring triangle.
            V_edge_X(i,j)=1/2*(V_X(triangles(i,mod(j,3)+1))+...
                V_X(triangles(i,mod(j+1,3)+1)));
            V_edge_Y(i,j)=1/2*(V_Y(triangles(i,mod(j,3)+1))+...
                V_Y(triangles(i,mod(j+1,3)+1)));
        else
            %Determine the location of the intersection of the triangle
            %edge and the line connecting the neighbouring triangle centers
            
            %Check if we already calculated the V_edge
            if (Tri_nb(i,j)<i)
                V_shared=triangles(i,:);
                V_shared(j)=[];
                %Filter out the vertices that are shared
                V_nb=triangles(Tri_nb(i,j),:);
                V_nb(V_nb==V_shared(1))=0;
                V_nb(V_nb==V_shared(2))=0;
                j2=find(V_nb);
                V_edge_X(i,j)=V_edge_X(Tri_nb(i,j),j2);
                V_edge_Y(i,j)=V_edge_Y(Tri_nb(i,j),j2);
            else 
                A=[V_X(triangles(i,mod(j+1,3)+1))-...
                    V_X(triangles(i,mod(j,3)+1)),...
                   -(V_c_X(Tri_nb(i,j))-V_c_X(i));...
                   V_Y(triangles(i,mod(j+1,3)+1))-...
                   V_Y(triangles(i,mod(j,3)+1)),...
                   -(V_c_Y(Tri_nb(i,j))-V_c_Y(i))];
                b=[V_c_X(i)-V_X(triangles(i,mod(j,3)+1));...
                   V_c_Y(i)-V_Y(triangles(i,mod(j,3)+1))];
                S=A\b;
                s=S(1);
                V_edge_X(i,j)=(1-s)*V_X(triangles(i,mod(j,3)+1))...
                             +s*V_X(triangles(i,mod(j+1,3)+1));
                V_edge_Y(i,j)=(1-s)*V_Y(triangles(i,mod(j,3)+1))...
                             +s*V_Y(triangles(i,mod(j+1,3)+1));  
            end
        end
    end 
end

if fig==1
    figure(1)
    hold on
    scatter(V_X,V_Y,'k')
    scatter(V_c_X,V_c_Y,'b')
    scatter(V_edge_X(:,1),V_edge_Y(:,1),'r')
    scatter(V_edge_X(:,2),V_edge_Y(:,2),'r')
    scatter(V_edge_X(:,3),V_edge_Y(:,3),'r')
    hold off
end
clear A b i j s S V_c

%% Construct the almost optimal PS triangles

%Initialise the vertices Q of the PS triangles
Q_X=zeros(Num_Ver,3);
Q_Y=zeros(Num_Ver,3);

%Function for determining if the centre of the polygon is in the generated
%triangle for the case of determining a triangle from 3 edges of the convex
%hull
signf  = @(p1,p2,p3)  (p1(1)- p3(1))*(p2(2)-p3(2))...
                     -(p2(1)- p3(1))*(p1(2)-p3(2));
%Source: 
%https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is
%-in-a-2d-triangle

if fig==1
    figure(1)
    hold on
end

for v=1:Num_Ver     %loop over triangles
    %Triangles that have v as a vertex
    [Tri_v,Ind_v]=find(triangles.ConnectivityList==v);
    
    %Preallocate PS-points
    PS_X=zeros(3*length(Tri_v)+1,1);
    PS_Y=zeros(3*length(Tri_v)+1,1);
    PS_X(1)=V_X(v);
    PS_Y(1)=V_Y(v);
    %Calculate PS-Points
    for i=1:length(Tri_v)
        PS_X(3*i-2+(1:3))=...
              [1/2*(V_edge_X(Tri_v(i),mod(Ind_v(i)-2,3)+1)+V_X(v));...
               1/2*(V_c_X(Tri_v(i))+V_X(v));...
               1/2*(V_edge_X(Tri_v(i),mod(Ind_v(i),3)+1)+V_X(v))];
        PS_Y(3*i-2+(1:3))=...
              [1/2*(V_edge_Y(Tri_v(i),mod(Ind_v(i)-2,3)+1)+V_Y(v));...
               1/2*(V_c_Y(Tri_v(i))+V_Y(v));...
               1/2*(V_edge_Y(Tri_v(i),mod(Ind_v(i),3)+1)+V_Y(v))];
    end
    
    %Filter the relevant points for the convex hull
    cvh=convhull(PS_X,PS_Y,'simplify',true);
    PS_X=PS_X(cvh);
    PS_Y=PS_Y(cvh);

    %Plot PS_points hull edges
%     scatter(PS_X,PS_Y)
    
    %Initialise the surface as infinity
    Sur=Inf;
    
    if fig==1
        figure(1)
        hold on
    end
    
    %loop over all triangles with 3 edges shared with the convex hull
    for i=1:length(cvh)-3
        u1=[PS_X(i);PS_Y(i)];
        u2=[PS_X(i+1);PS_Y(i+1)];
        for j=i+1:length(cvh)-2
            v1=[PS_X(j);PS_Y(j)];
            v2=[PS_X(j+1);PS_Y(j+1)]; 

            A1=[u2-u1,v1-v2];
            if (abs(det(A1))>10^-14)
                S=A1\(v1-u1);
                Q_temp(:,1)=u1+S(1)*(u2-u1);
                for k=j+1:length(cvh)-1
                    w1=[PS_X(k);PS_Y(k)];
                    w2=[PS_X(k+1);PS_Y(k+1)];

                    A2=[v2-v1,w1-w2];
                    if (abs(det(A2))>10^-14)
                        S=A2\(w1-v1);
                        Q_temp(:,2)=v1+S(1)*(v2-v1);

                        A3=[w2-w1,u1-u2];
                        if (abs(det(A3))>10^-14)
                            S=A3\(u1-w1);
                            Q_temp(:,3)=w1+S(1)*(w2-w1);

                            %Prepare triangle validness check
                            b1 = signf([V_X(v);V_Y(v)],...
                                Q_temp(:,1), Q_temp(:,2)) < 0;
                            b2 = signf([V_X(v);V_Y(v)],...
                                Q_temp(:,2), Q_temp(:,3)) < 0;
                            b3 = signf([V_X(v);V_Y(v)],...
                                Q_temp(:,3), Q_temp(:,1)) < 0;

                            %Check if center point is in the triangle
                            if ((b1 == b2) && (b2 == b3))
                                %Check if the surface is smaller than the
                                %previous best
                                [~,Sur_temp]=convhull(Q_temp');
                                if (Sur_temp<Sur)
                                    Sur=Sur_temp;
                                    Q_X(v,:)=Q_temp(1,:);
                                    Q_Y(v,:)=Q_temp(2,:);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %loop over all triangles with 2 edges shared with the convex hull
    hold on
    for i=1:length(cvh)-2
        u1=[PS_X(i);PS_Y(i)];
        u2=[PS_X(i+1);PS_Y(i+1)];
        for j=i+1:length(cvh)-1
            v1=[PS_X(j);PS_Y(j)];
            v2=[PS_X(j+1);PS_Y(j+1)]; 

            A1=[u2-u1,v1-v2];
            if (abs(det(A1))>10^-13)
                %First, find the first SP-triangle vertex Q1
                S=A1\(v1-u1);
                Q_temp(:,1)=u1+S(1)*(u2-u1);

                %Then find the bisector corresponding to that PS-triangle
                %vertex
                u_norm=(1/2*(u2+u1)-Q_temp(:,1))./...
                    norm(1/2*(u2+u1)-Q_temp(:,1));
                v_norm=(1/2*(v2+v1)-Q_temp(:,1))./...
                    norm(1/2*(v2+v1)-Q_temp(:,1));
                UV_bis=u_norm+v_norm;

                %Determine which PS point is furthest away from Q1 in the
                %inner product semi-norm with the bisector of the angle of
                %Q1
                PS_Q_vectors=[PS_X(1:end-1)';PS_Y(1:end-1)']...
                                -Q_temp(:,1)*ones(1,length(cvh)-1);
                Distances=dot(PS_Q_vectors,UV_bis*ones(1,length(cvh)-1));
                PS_max_dist_ind=find(Distances==max(Distances),1);

                PS_max_dist=[PS_X(PS_max_dist_ind);PS_Y(PS_max_dist_ind)]; 

                %Now determine the vector orthogonal to the bisector
                %through Q1
                UV_ort=[-UV_bis(2);UV_bis(1)];

                %Finally, determine the Q2 and Q3 by intersecting the line
                %constructed with UV_ort through NS_max_dist with the other
                %two lines through Q1
                A2=[u2-u1,UV_ort];
                S=A2\(PS_max_dist-u1);
                Q_temp(:,2)=u1+S(1)*(u2-u1);

                A3=[v2-v1,UV_ort];
                S=A3\(PS_max_dist-v1);
                Q_temp(:,3)=v1+S(1)*(v2-v1);

                %Check if the surface is smaller than previous best
                [~,Sur_temp]=convhull(Q_temp');
                if (Sur_temp<Sur)
                    Sur=Sur_temp;
                    Q_X(v,:)=Q_temp(1,:);
                    Q_Y(v,:)=Q_temp(2,:);
                end
            end
        end
    end

%     %Uncomment for the location of the PS-triangle vertices
    if fig==1
        scatter(Q_X(v,:),Q_Y(v,:),'filled')
    end
    
end
hold off

clear A1 A2 A3 b1 b2 b3 cvh Distances i Ind_v j k PS_max_dist ...
    PS_max_dist_ind PS_Q_vectors PS_X PS_Y Q_temp S signf Sur Sur_temp ...
    Tri_v u1 u2 u_norm UV_bis UV_ort v v1 v2 v_norm w1 w2 
clear fig

%% Calculating Powell-Sabin Spline coefficients

%This part is based on the paper of Paul Dierckx, 
%"On calculating normalized Powell-Sabin B-splines"

%First calculate the alphas, betas and gammas corresponding to the
%coordinates of the BS-triangle (Equations (51) and (52) from Dierckx

%alpha are the barycentric coordinates of V1 w.r.t. the PS triangle
alpha=zeros(Num_Ver,3);
e    =zeros(Num_Ver,1);
beta =zeros(Num_Ver,3);
gamma=zeros(Num_Ver,3);

for i=1:Num_Ver
    %Equation 45
    alpha(i,:)=[Q_X(i,:);Q_Y(i,:);[1 1 1]]\[V_X(i);V_Y(i);1];
    
    %Equation 44        
    e(i)=det([Q_X(i,:)',Q_Y(i,:)',[1; 1; 1]]);
    %Equation 51,52
    beta(i,:) =1/e(i)*...
        [Q_Y(i,2)-Q_Y(i,3),Q_Y(i,3)-Q_Y(i,1),Q_Y(i,1)-Q_Y(i,2)];
    gamma(i,:)=1/e(i)*...
        [Q_X(i,3)-Q_X(i,2),Q_X(i,1)-Q_X(i,3),Q_X(i,2)-Q_X(i,1)];  
end

clear i

%% Define the basis functions by their Bezier Ordinates on all subtriangles

%The figures and equations refered to can be found in  
%"On calculating normalized Powell-Sabin B-splines", by Paul Dierckx.

%Preallocate coordinates of the Bezier Ordinates
B_o_X=zeros(Num_Tri,6,6);
B_o_Y=zeros(Num_Tri,6,6);

%X and Y locations corresponding Bezier Ordinates for each main triangle
for tr=1:Num_Tri
    v1=triangles(tr,1);
    v2=triangles(tr,2);
    v3=triangles(tr,3);
    
    B_o_X(tr,:,:)=[V_c_X(tr) 1/2*(V_c_X(tr)+V_X(v1)) V_X(v1) ...
                       1/2*(V_X(v1)+V_edge_X(tr,3)) V_edge_X(tr,3) ...
                       1/2*(V_edge_X(tr,3)+V_c_X(tr));...
                   V_c_X(tr) 1/2*(V_c_X(tr)+V_edge_X(tr,3)) ...
                       V_edge_X(tr,3) 1/2*(V_edge_X(tr,3)+V_X(v2)) ...
                       V_X(v2) 1/2*(V_X(v2)+V_c_X(tr));...
                   V_c_X(tr) 1/2*(V_c_X(tr)+V_X(v2)) V_X(v2) ...
                       1/2*(V_X(v2)+V_edge_X(tr,1)) V_edge_X(tr,1) ...
                       1/2*(V_edge_X(tr,1)+V_c_X(tr));...
                   V_c_X(tr) 1/2*(V_c_X(tr)+V_edge_X(tr,1)) ...
                       V_edge_X(tr,1) 1/2*(V_edge_X(tr,1)+V_X(v3)) ...
                       V_X(v3) 1/2*(V_X(v3)+V_c_X(tr));...
                   V_c_X(tr) 1/2*(V_c_X(tr)+V_X(v3)) V_X(v3) ...
                       1/2*(V_X(v3)+V_edge_X(tr,2)) V_edge_X(tr,2) ...
                       1/2*(V_edge_X(tr,2)+V_c_X(tr));...
                   V_c_X(tr) 1/2*(V_c_X(tr)+V_edge_X(tr,2)) ...
                       V_edge_X(tr,2) 1/2*(V_edge_X(tr,2)+V_X(v1))...
                       V_X(v1) 1/2*(V_X(v1)+V_c_X(tr))];
    
    B_o_Y(tr,:,:)=[V_c_Y(tr) 1/2*(V_c_Y(tr)+V_Y(v1)) V_Y(v1) ...
                       1/2*(V_Y(v1)+V_edge_Y(tr,3)) V_edge_Y(tr,3) ...
                       1/2*(V_edge_Y(tr,3)+V_c_Y(tr));...
                   V_c_Y(tr) 1/2*(V_c_Y(tr)+V_edge_Y(tr,3))  ...
                       V_edge_Y(tr,3) 1/2*(V_edge_Y(tr,3)+V_Y(v2)) ...
                       V_Y(v2) 1/2*(V_Y(v2)+V_c_Y(tr));...
                   V_c_Y(tr) 1/2*(V_c_Y(tr)+V_Y(v2)) V_Y(v2) ...
                       1/2*(V_Y(v2)+V_edge_Y(tr,1)) V_edge_Y(tr,1) ...
                       1/2*(V_edge_Y(tr,1)+V_c_Y(tr));...
                   V_c_Y(tr) 1/2*(V_c_Y(tr)+V_edge_Y(tr,1)) ...
                       V_edge_Y(tr,1) 1/2*(V_edge_Y(tr,1)+V_Y(v3)) ...
                       V_Y(v3) 1/2*(V_Y(v3)+V_c_Y(tr));...
                   V_c_Y(tr) 1/2*(V_c_Y(tr)+V_Y(v3)) V_Y(v3) ...
                       1/2*(V_Y(v3)+V_edge_Y(tr,2)) V_edge_Y(tr,2) ...
                       1/2*(V_edge_Y(tr,2)+V_c_Y(tr));...
                   V_c_Y(tr) 1/2*(V_c_Y(tr)+V_edge_Y(tr,2))  ...
                       V_edge_Y(tr,2) 1/2*(V_edge_Y(tr,2)+V_Y(v1)) ...
                       V_Y(v1) 1/2*(V_Y(v1)+V_c_Y(tr))];
end
            
%Preallocate Bezier Ordinates. B_O(Ver,Q,Tri,Subtri,BO) denotes the basis
%function belonging to PS triangle vertex Q of main vertex Ver. Then for
%each triangle in the mesh (1 to Num_Tri), each of its subtriangles Sub_Tri
%(1 to 6) is considered, and of each subtriangle, each of the Bezier-
%Ordinates BO (1 to 6) is given the appropriate value to the basis
%function. Most of this 5D matrix will be empty due to the finite support
%of each basis function: only the triangles that have Ver as one of their
%vertices will have nonzero entries. Unfortunately, Matlab does not allow
%standard to convert a 5D matrix to sparse.
%B_O(Ver,Q,Tri,Subtri,BO)
B_o  =zeros(Num_Ver,3,Num_Tri,6,6);

for v=1:Num_Ver
    for q=1:3
        %Find the triangles with vertex v and also save the index of vertex
        %v in this triangle (1, 2 or 3).
        [tri_v,ind_v]=find(triangles.ConnectivityList==v);
        for t=1:length(tri_v)
            %Number and index of the current main triangle
            tr=tri_v(t);
            ind=ind_v(t);
            
            %Barycentric coordinates of the center vertex
            BC_c=cartesianToBarycentric(triangles,tr,...
                [V_c_X(tr),V_c_Y(tr)]);
            %Keep in mind that we are treating vertex ind of Tri as number
            %V1 in Fig.3 here. E.g.: Tri=(5,4,6), and we are looking at
            %vertex 4, then ind=2, then (V1,V2,V3) of Fig. 3 is
            %(V4,V6,V5) in this example. 
            a=BC_c(ind);
            b=BC_c(mod(ind,3)+1);
            c=BC_c(mod(ind+1,3)+1);
            
            %Barycentric coordinates of the edge vertices
            BC_edge=cartesianToBarycentric(triangles,[tr;tr;tr],...
                [V_edge_X(tr,:)',V_edge_Y(tr,:)']); 
            
            %Barycentric coordinates of the edge vertices, some of these
            %are obsolete. More information can be found at the bottom of
            %page 65,
            lambda1=BC_edge(mod(ind+1,3)+1,ind);
%             lambda2=BC_edge(mod(ind+1,3)+1,mod(ind,3)+1);
%             mu2    =BC_edge(ind           ,mod(ind,3)+1);
%             mu3    =BC_edge(ind           ,mod(ind+1,3)+1);
            nu1    =BC_edge(mod(ind,3)+1  ,ind);
%             nu3    =BC_edge(mod(ind,3)+1  ,mod(ind+1,3)+1);
            
            %The locations of the three vertices of triangle Tri, counted
            %counterclockwise, with vertex v the first one. 
            %Please see Fig. 3
            x1=V_X(v);
            x2=V_X(triangles.ConnectivityList(tr,mod(ind,3)+1));
            x3=V_X(triangles.ConnectivityList(tr,mod(ind+1,3)+1));
            y1=V_Y(v);
            y2=V_Y(triangles.ConnectivityList(tr,mod(ind,3)+1));
            y3=V_Y(triangles.ConnectivityList(tr,mod(ind+1,3)+1));
            
            %Equation 15,16
            beta_bar =beta(v,q)*(x2-x1)+gamma(v,q)*(y2-y1);
            gamma_bar=beta(v,q)*(x3-x1)+gamma(v,q)*(y3-y1);
            
            %Equation 12,13,14
            L=alpha(v,q)+(1-lambda1)/2*beta_bar;
            L_prime=alpha(v,q)+(1-nu1)/2*gamma_bar;
            L_tilde=alpha(v,q)+b/2*beta_bar+c/2*gamma_bar;

            %Define the Bezier ordinates (Fig. 4)
            B_O_temp=[a*L_tilde L_tilde alpha(v,q) L ...
                        lambda1*L lambda1*L_tilde;...
                      a*L_tilde lambda1*L_tilde lambda1*L 0 0 0;...
                      a*L_tilde 0 0 0 0 0; ...
                      a*L_tilde 0 0 0 0 0; ...
                      a*L_tilde 0 0 0 nu1*L_prime nu1*L_tilde;...
                      a*L_tilde nu1*L_tilde nu1*L_prime ... 
                        L_prime alpha(v,q) L_tilde];
            
            %Reorder B_O_temp such that it fits with the right vertex.
            B_O_temp=circshift(B_O_temp,2*(ind-1),1);
            
            %Insert B_O_temp on the right position.
            B_o(v,q,tr,:,:)=B_O_temp;
        end
    end
end

B_o(B_o<0)=0;
clear a b B_O_temp BC_c BC_edge beta_bar c gamma_bar ind ind_v L ...
    L_prime L_tilde lambda1 lambda2 mu2 mu3 nu1 nu3 q t tr tri_v v ...
    x1 x2 x3 y1 y2 y3 Num_tri
     
end
