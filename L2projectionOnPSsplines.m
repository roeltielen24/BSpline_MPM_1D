%This code has been designed to find a L2-projection of an arbitrary 
%function. The function will be projected on a base of piecewise quadratic
%B-splines with local support. This code was written to test the base
%functions, while not having to worry about boundary conditions yet. Thus
%the choice for L2-projection instead of the Poisson function or any more
%difficult (P)DE.
%The function is approximated by solving
%
%\int phi_i(f*-f)dOmega=0, for all phi_j; f*=sum_j c_j*phi_j
%
%This leads to the system of equations
%
%sum_j c_i[int(phi_i*phi_j dOmega)] = int(phi_i*f dOmega), for i=1...m 
%

clc
clear
% close all

%% Parameters (This part, and only this part, is for tweaking!)
%Computation time measurement
tic
disp('Initialising...')

%Dimension of the problem, do not change this, the code is specialised for
%2D!
dim=2;

%Enable lumping with 1, no lumping with 0, and 3x3 block lumping with 2
lumping = 2;

%Function to be approximated
% f=@(x,y) 0.1*sin(pi.*x);
f=@(x,y) 1.*y;
% f=@(x,y) heaviside(x-0.7).*(x-0.7).*heaviside(y-0.7).*(y-0.7);

%In case we want to impose a essential homogeneous boundary condition on
%all sides in correspondence with the approximated function, put BC
%(Boundary Condition) to 1.
BC=0;

%Number of grid points in each direction (at least 2). 
n_x=4+1;
n_y=4+1;
%!!!Note: there will be (3*n) basis functions, so the Mass Matrix will be
%(3*n)^2 times (3*n)^2. Therefore, it is advised to keep n smaller than 20.

%Choose for a regular (=1) or randomly generated (=0) grid
regular=1;

%Gauss integration degree, choose 1, 2, 3, 4 or 5 (6 and higher not
%implemented). Gauss degree 4 yields exact integration of the Mass matrix.
Gauss_degree=5;

%Resolution of the subtriangles. Note that there are about 1.5*n^2 main
%triangles, and each main triangle has 6 subtriangles. Therefore, the
%number of data points will be about 9*n^2*resolution
resolution=10;  %Advised range: 3-20. Reduce resolution for large n. 
resolution=resolution+1;

toc
%% Setting up grid
tic
disp('Setting up...')

if regular==1
    %Perfecly divided in squares, V for Vertex
    V_X=linspace(0,1,n_x)'*ones(1,n_y);
    V_X=reshape(V_X,[n_x*n_y,1]);
    V_Y=ones(n_x,1)*linspace(0,1,n_y);
    V_Y=reshape(V_Y,[n_x*n_y,1]);
    %Vertex Numbering goes as follows:
    %13 14 15 16
    % 9 10 11 12
    % 5  6  7  8
    % 1  2  3  4
else
    %Irregular polygon
    %The center part is randomly generated, the edges are regular
    V_X=[zeros(1,n_y-2);rand(n_x-2,n_y-2);ones(1,n_y-2)];
    V_X=[linspace(0,1,n_x)',V_X,linspace(0,1,n_x)'];
    V_X=reshape(V_X,[n_x*n_y,1]);
    V_Y=[zeros(n_x-2,1),rand(n_x-2,n_y-2),ones(n_x-2,1)];
    V_Y=[linspace(0,1,n_y);V_Y;linspace(0,1,n_y)];
    V_Y=reshape(V_Y,[n_x*n_y,1]);
end

%Number of vertices
Num_Ver=length(V_X);

%Generate triangles
Tri = delaunayTriangulation(V_X,V_Y);
Num_Tri=size(Tri,1);

%Generate a list of neighours
Tri_nb = neighbors(Tri);

%Plot the mesh
figure(1)
hold on
trimesh(Tri,V_X,V_Y)
hold off
title('Triangulations')

toc
%% Construct Powell Sabin refinement
tic
disp('Constructing PS refinement...')

%Generate points in the center of each triangle with incenter, which is the
%intersection point of the three bisectors of the angles of the triangle. 
%http://mathworld.wolfram.com/Incenter.html
V_c=incenter(Tri);
V_c_X=V_c(:,1);
V_c_Y=V_c(:,2);

figure(1)
hold on
% scatter(V_c(:,1),V_c(:,2))
hold off

%Define V_edge_X/Y(Tri,Ver) as the vertex in Triangle Tri on the OPPOSITE
%side of the Vertex it belongs to. This way, the numbering is also
%unambiguous for tetrahedrons.
V_edge_X=zeros(Num_Tri,3);
V_edge_Y=zeros(Num_Tri,3);

%Determine the inner triangle points on the edges of the main triangle.
for i=1:Num_Tri
    for j=1:3   %Number of edges
        if isnan(Tri_nb(i,j)) %No neighbouring triangle.
            V_edge_X(i,j)=...
                1/2*(V_X(Tri(i,mod(j,3)+1))+V_X(Tri(i,mod(j+1,3)+1)));
            V_edge_Y(i,j)=...
                1/2*(V_Y(Tri(i,mod(j,3)+1))+V_Y(Tri(i,mod(j+1,3)+1)));
        else
            %Determine the location of the intersection of the triangle
            %edge and the line connecting the neighbouring triangle centers
            A=[V_X(Tri(i,mod(j+1,3)+1))-V_X(Tri(i,mod(j,3)+1)),...
               -(V_c_X(Tri_nb(i,j))-V_c_X(i));...
               V_Y(Tri(i,mod(j+1,3)+1))-V_Y(Tri(i,mod(j,3)+1)),...
               -(V_c_Y(Tri_nb(i,j))-V_c_Y(i))];
            b=[V_c_X(i)-V_X(Tri(i,mod(j,3)+1));...
               V_c_Y(i)-V_Y(Tri(i,mod(j,3)+1))];
            S=A\b;
            s=S(1);
            V_edge_X(i,j)=(1-s)*V_X(Tri(i,mod(j,3)+1))...
                         +s*V_X(Tri(i,mod(j+1,3)+1));
            V_edge_Y(i,j)=(1-s)*V_Y(Tri(i,mod(j,3)+1))...
                         +s*V_Y(Tri(i,mod(j+1,3)+1));  
        end
    end 
end

figure(1)
hold on
xlabel('x')
ylabel('y')
scatter(V_X,V_Y,'k')
scatter(V_c_X,V_c_Y,'b')
scatter(V_edge_X(:,1),V_edge_Y(:,1),'r')
scatter(V_edge_X(:,2),V_edge_Y(:,2),'r')
scatter(V_edge_X(:,3),V_edge_Y(:,3),'r')
hold off

clear A b i j s S V_c
toc
%% Construct the almost optimal PS triangles
tic
disp('Calculating almost optimal PS triangles...')

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

figure(1)
hold on
for v=1:Num_Ver     %loop over triangles
    %Triangles that have v as a vertex
    [Tri_v,Ind_v]=find(Tri.ConnectivityList==v);
    
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

    figure(1)
    hold on
    %loop over all triangles with 3 edges shared with the convex hull
    for i=1:length(cvh)-3
        u1=[PS_X(i);PS_Y(i)];
        u2=[PS_X(i+1);PS_Y(i+1)];
        for j=i+1:length(cvh)-2
            v1=[PS_X(j);PS_Y(j)];
            v2=[PS_X(j+1);PS_Y(j+1)]; 

            A1=[u2-u1,v1-v2];
            if (abs(det(A1))>10^-15)
                S=A1\(v1-u1);
                Q_temp(:,1)=u1+S(1)*(u2-u1);
                for k=j+1:length(cvh)-1
                    w1=[PS_X(k);PS_Y(k)];
                    w2=[PS_X(k+1);PS_Y(k+1)];

                    A2=[v2-v1,w1-w2];
                    if (abs(det(A2))>10^-15)
                        S=A2\(w1-v1);
                        Q_temp(:,2)=v1+S(1)*(v2-v1);

                        A3=[w2-w1,u1-u2];
                        if (abs(det(A3))>10^-15)
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
            if (abs(det(A1))>10^-14)
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
    
    scatter(Q_X(v,:),Q_Y(v,:),'filled')
    
end
hold off

clear A1 A2 A3 b1 b2 b3 cvh Distances i Ind_v j k PS_max_dist ...
    PS_max_dist_ind PS_Q_vectors PS_X PS_Y Q_temp S signf Sur Sur_temp ...
    Tri_v u1 u2 u_norm UV_bis UV_ort v v1 v2 v_norm w1 w2 
toc
%% Calculating Powell-Sabin Spline coefficients
tic
disp('Calculating PS spline coefficients...')

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
toc
%% Define the basis functions by their Bezier Ordinates on all subtriangles
tic
disp('Calculating Bezier Ordinates of the basis functions...')

%The figures and equations refered to can be found in  
%"On calculating normalized Powell-Sabin B-splines", by Paul Dierckx.

%Preallocate coordinates of the Bezier Ordinates
B_O_X=zeros(Num_Tri,6,6);
B_O_Y=zeros(Num_Tri,6,6);

%X and Y locations corresponding Bezier Ordinates for each main triangle
for tr=1:Num_Tri
    v1=Tri(tr,1);
    v2=Tri(tr,2);
    v3=Tri(tr,3);
    
    B_O_X(tr,:,:)=[V_c_X(tr) 1/2*(V_c_X(tr)+V_X(v1)) V_X(v1) ...
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
    
    B_O_Y(tr,:,:)=[V_c_Y(tr) 1/2*(V_c_Y(tr)+V_Y(v1)) V_Y(v1) ...
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
B_O  =zeros(Num_Ver,3,Num_Tri,6,6);

for v=1:Num_Ver
    for q=1:3
        %Find the triangles with vertex v and also save the index of vertex
        %v in this triangle (1, 2 or 3).
        [tri_v,ind_v]=find(Tri.ConnectivityList==v);
        for t=1:length(tri_v)
            %Number and index of the current main triangle
            tr=tri_v(t);
            ind=ind_v(t);
            
            %Barycentric coordinates of the center vertex
            BC_c=cartesianToBarycentric(Tri,tr,[V_c_X(tr),V_c_Y(tr)]);
            %Keep in mind that we are treating vertex ind of Tri as number
            %V1 in Fig.3 here. E.g.: Tri=(5,4,6), and we are looking at
            %vertex 4, then ind=2, then (V1,V2,V3) of Fig. 3 is
            %(V4,V6,V5) in this example. 
            a=BC_c(ind);
            b=BC_c(mod(ind,3)+1);
            c=BC_c(mod(ind+1,3)+1);
            
            %Barycentric coordinates of the edge vertices
            BC_edge=cartesianToBarycentric(Tri,[tr;tr;tr],...
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
            x2=V_X(Tri.ConnectivityList(tr,mod(ind,3)+1));
            x3=V_X(Tri.ConnectivityList(tr,mod(ind+1,3)+1));
            y1=V_Y(v);
            y2=V_Y(Tri.ConnectivityList(tr,mod(ind,3)+1));
            y3=V_Y(Tri.ConnectivityList(tr,mod(ind+1,3)+1));
            
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
            B_O(v,q,tr,:,:)=B_O_temp;
        end
    end
end

clear a b B_O_temp BC_c BC_edge beta_bar c gamma_bar ind ind_v L ...
    L_prime L_tilde lambda1 lambda2 mu2 mu3 nu1 nu3 q t tr tri_v v ...
    x1 x2 x3 y1 y2 y3
toc
%% Construct the Mass matrix M and the Right Hand Side (RHS)
tic
disp('Constructing Mass matrix and Right Hand Side...')
%The figures and equations refered to can be found in  
%"On calculating normalized Powell-Sabin B-splines", by Paul Dierckx.

%Preallocate the mass matrix. M_ij contains \int phi_i*phi_j dOmega, where
%i=3*(v-1)+q, with v the vertex number (from 1 to Num_Ver) and q the number
%of the vertex of the PS triangle (from 1 to 3).
M=zeros(3*Num_Ver);

%Preallocate the RHS
RHS=zeros(3*Num_Ver,1);

%Gau\ss quadrature rule for triangles in barycentric coordinates. Source:
%http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
if Gauss_degree==1
    %Degree 1
    Gauss_t1=1/3;
    Gauss_t2=1/3;
    Gauss_t3=1/3;
    Gauss_W =1;
elseif Gauss_degree==2
    %Degree 2
    Gauss_t1=[2/3,1/6,1/6];
    Gauss_t2=[1/6,2/3,1/6];
    Gauss_t3=[1/6,1/6,2/3];
    Gauss_W =[1/3,1/3,1/3];
elseif Gauss_degree==3
    %Degree 3
    Gauss_t1=     [ 1/3,3/5,1/5,1/5];
    Gauss_t2=     [ 1/3,1/5,3/5,1/5];
    Gauss_t3=     [ 1/3,1/5,1/5,3/5];
    Gauss_W =1/48*[-27 ,25 ,25 ,25 ];
elseif Gauss_degree==4
    %Degree 4
    Gauss_t1=[0.445948490915970 0.445948490915970 0.108103018168070...
        0.0915762135097700 0.0915762135097700 0.816847572980460];
    Gauss_t2=[0.445948490915970 0.108103018168070 0.445948490915970...
        0.0915762135097700 0.816847572980460 0.0915762135097700];
    Gauss_t3=1-Gauss_t1-Gauss_t2;
    Gauss_W =[0.223381589678010 0.223381589678010 0.223381589678010...
        0.109951743655320 0.109951743655320 0.109951743655320];
elseif Gauss_degree==5
    %Degree 5
    Gauss_t1=[0.333333333333330  0.470142064105110 0.470142064105110...
              0.0597158717897700 0.101286507323460 0.101286507323460...
              0.797426985353090];
    Gauss_t2=[0.333333333333330 0.470142064105110 0.0597158717897700...
              0.470142064105110 0.101286507323460 0.797426985353090...
              0.101286507323460];
    Gauss_t3=1-Gauss_t1-Gauss_t2;
    Gauss_W =[0.225000000000000 0.132394152788510 0.132394152788510... 
              0.132394152788510 0.125939180544830 0.125939180544830... 
              0.125939180544830];
else 
    disp('Chosen Gauss degree is not implemented,')
    disp('choose 1, 2, 3, 4 or 5 instead.')
    return
end

%For Gau\ss integrating each subtriangle, the surface of each subtriangle
%is required. Here St_sur(tr,s) represents the surface of the subtriangle s
%(from 1 to 6) of main triangle s (from 1 to Num_Tri).
%St_sur(Tri,Subtri)
St_sur=zeros(Num_Tri,6);

%loop over all i (rows) of M_ij and RHS_i
%loop over all vertices
for vi=1:Num_Ver
    %loop over all vertices of the PS triangle corresponing to v
    for qi=1:3
        %i for M_ij and RHS_i.
        i=3*(vi-1)+qi;
        
        %loop over all j (columns) of M_ij
        %We only start at vertex vi, because M is symmetric, so we can copy
        %the lower triangle from the upper one
        for vj=vi:Num_Ver
            for qj=1:3
                %j for M_ij.
                j=3*(vj-1)+qj;
                
                %Find all triangles that have both vi and vj as vertices
                Tri_vi=sum(Tri.ConnectivityList==vi,2);
                Tri_vj=sum(Tri.ConnectivityList==vj,2);
                Tri_vivj=find(Tri_vi==1 & Tri_vj==1);  
                
                M_temp=0;
                %loop over all triangles that have both vi and vj as 
                %vertices
                for tr=Tri_vivj'
                    
                    %loop over all subtriangles in each triangle
                    for s=1:6
                        %The surface of the triangle, necessary to do Gauss
                        %integration
                        [~,St_sur(tr,s)]=convhull(...
                            squeeze(B_O_X(tr,s,:)),squeeze(B_O_Y(tr,s,:)));   
                        
                        %loop over all Gauss points
                        for gp=1:length(Gauss_W)
                            %Gauss_weight*(phi_i(gp)*phi_j(gp))
                            M_temp=M_temp+Gauss_W(gp)*St_sur(tr,s)*(...
                                B_O(vi,qi,tr,s,1)*Gauss_t1(gp)^2+...
                                B_O(vi,qi,tr,s,2)*...
                                    2*Gauss_t1(gp)*Gauss_t2(gp)+...
                                B_O(vi,qi,tr,s,3)*Gauss_t2(gp)^2+...
                                B_O(vi,qi,tr,s,4)*...
                                    2*Gauss_t2(gp)*Gauss_t3(gp)+...
                                B_O(vi,qi,tr,s,5)*Gauss_t3(gp)^2+...
                                B_O(vi,qi,tr,s,6)*...
                                    2*Gauss_t3(gp)*Gauss_t1(gp) )...
                                * (...
                                B_O(vj,qj,tr,s,1)*Gauss_t1(gp)^2+...
                                B_O(vj,qj,tr,s,2)*...
                                    2*Gauss_t1(gp)*Gauss_t2(gp)+...
                                B_O(vj,qj,tr,s,3)*Gauss_t2(gp)^2+...
                                B_O(vj,qj,tr,s,4)*...
                                    2*Gauss_t2(gp)*Gauss_t3(gp)+...
                                B_O(vj,qj,tr,s,5)*Gauss_t3(gp)^2+...
                                B_O(vj,qj,tr,s,6)*...
                                    2*Gauss_t3(gp)*Gauss_t1(gp) );
                        end
                    end
                    
                end
                
                %Insert M_temp into M
                M(i,j)=M_temp;
                %Due to symmetry M_ij=M_ji
                M(j,i)=M_temp;
            end
        end
        
        %Calculate RHS_temp
        RHS_temp=0;
        
        %Find all triangles that have vi as vertices
        Tri_vi=sum(Tri.ConnectivityList==vi,2);
        Tri_vi=find(Tri_vi==1);
        
        %loop over all triangles that have both vi as vertices
        for tr=Tri_vi'
        
            %loop over all subtriangles in each triangle
            for s=1:6
                %The surface of the triangle, necessary to do Gauss
                %integration
                [~,St_sur(tr,s)]=convhull(...
                    squeeze(B_O_X(tr,s,:)),squeeze(B_O_Y(tr,s,:)));   

                %loop over all Gauss points
                for gp=1:length(Gauss_W)
                    %X and Y location of this Gauss point
                    gp_X=Gauss_t1(gp)*B_O_X(tr,s,1)+...
                         Gauss_t2(gp)*B_O_X(tr,s,3)+...
                         Gauss_t3(gp)*B_O_X(tr,s,5);
                    gp_Y=Gauss_t1(gp)*B_O_Y(tr,s,1)+...
                         Gauss_t2(gp)*B_O_Y(tr,s,3)+...
                         Gauss_t3(gp)*B_O_Y(tr,s,5);                     
                    %Gauss_weight*(phi_i(gp)*phi_j(gp))
                    RHS_temp=RHS_temp+Gauss_W(gp)*St_sur(tr,s)*( ...
                        B_O(vi,qi,tr,s,1)*Gauss_t1(gp)^2+...
                        B_O(vi,qi,tr,s,2)*2*Gauss_t1(gp)*Gauss_t2(gp)+...
                        B_O(vi,qi,tr,s,3)*Gauss_t2(gp)^2+...
                        B_O(vi,qi,tr,s,4)*2*Gauss_t2(gp)*Gauss_t3(gp)+...
                        B_O(vi,qi,tr,s,5)*Gauss_t3(gp)^2+...
                        B_O(vi,qi,tr,s,6)*2*Gauss_t3(gp)*Gauss_t1(gp) )...
                        *...
                        ( f(gp_X,gp_Y) );
                end
            end

        end
        
        %Insert RHS_remp in RHS
        RHS(i)=RHS_temp;        
    end
end

clear gp gp_X gp_Y i j M_temp qi qj s tr Tri_vi Tri_vivj Tri_vj ...
    v1 v2 v3 vi vj    
toc


%% Implement boundary conditions

%In case of a homogeneous Dirichlet boundary condition on all sides
if BC==1
    Bound_south=find(abs(reshape(Q_Y',[],1)  )<1e-14);
    Bound_north=find(abs(reshape(Q_Y',[],1)-1)<1e-14);
    Bound_west= find(abs(reshape(Q_X',[],1)  )<1e-14);
    Bound_east= find(abs(reshape(Q_X',[],1)-1)<1e-14);
    
    Bound=union(union(union(Bound_south,Bound_north)...
        ,Bound_east),Bound_west);
    M(Bound,:)=0;
    M(:,Bound)=0;
    M(Bound,Bound)=eye(length(Bound));
    RHS(Bound_south)=0;
    RHS(Bound_north)=0;
    RHS(Bound_west)=0;
    RHS(Bound_east)=0;
end
if BC == 2
    Bound_west= find(abs(reshape(Q_X',[],1)  )<1e-14);
    Bound_east= find(abs(reshape(Q_X',[],1)-1)<1e-14);
    
    Bound=union(Bound_east,Bound_west);
    M(Bound,:)=0;
    M(:,Bound)=0;
    M(Bound,Bound)=eye(length(Bound));
    RHS(Bound_west)=0;
    RHS(Bound_east)=0;
end

%% Solve for the coefficients of the basis functions: Ac=b
tic
disp('Solving for the coefficients of the basis functions...')

if lumping == 1
    M_lumped = sum(M,2);
    c=RHS./M_lumped;
elseif lumping == 2
    M_blockLumped = zeros(size(M));
    M_copy=M;
    for i=1:n_x*n_y
        M_blockLumped(3*(i-1)+(1:3),3*(i-1)+(1:3)) = ...
            M(3*(i-1)+(1:3),3*(i-1)+(1:3));
        M_copy(3*(i-1)+(1:3),3*(i-1)+(1:3))=0;
    end
    M_blockLumped=M_blockLumped+diag(sum(M_copy,2));
    c=M_blockLumped\RHS;
else
    c=M\RHS;
end

% c=zeros(n_x*n_y,3); c(8,2)=1;

%Shape c such that c(Ver,Q) represents the coefficient of the basis 
%function of Vertex Ver and PS-triangle vertex Q.
c=reshape(c',[3,Num_Ver])';

toc
%% Plot the solution
tic
disp('Visualising the solution...')

%Generate Barycentric coordinates for the grid points in each subtriangle
BC_1=zeros(1,resolution*(resolution+1)/2);
BC_2=zeros(1,resolution*(resolution+1)/2);
%Counter
k=0;
for i=1:resolution
    for j=1:resolution-i+1
        k=k+1;
        BC_1(k)=(i-1)/(resolution-1);
        BC_2(k)=(j-1)/(resolution-1);
    end
end
BC_3=1-BC_1-BC_2;

%The X and Y location of the approximation Phi and the value of Phi, stored
%as follows. Phi(Tri,Subtri,loc) represents the value of Phi in the
%datapoint at triangle Tri, at its subtriangle Subtri, and more
%specifically at the predetermined barycentric coordinates loc
Phi_X=zeros(Num_Tri,6,resolution*(resolution+1)/2);
Phi_Y=zeros(Num_Tri,6,resolution*(resolution+1)/2);
Phi  =zeros(Num_Tri,6,resolution*(resolution+1)/2);

%Loop over all main triangles
for tr=1:Num_Tri
    
    %Loop over all subtriangles
    for s=1:6
        
        %loop over all locations in the subtriangle
        for loc=1:length(BC_1)
            
            %Determine the X and Y values of the current location
            Phi_X(tr,s,loc)=BC_1(loc)*B_O_X(tr,s,1)+...
                            BC_2(loc)*B_O_X(tr,s,3)+...
                            BC_3(loc)*B_O_X(tr,s,5);
            Phi_Y(tr,s,loc)=BC_1(loc)*B_O_Y(tr,s,1)+...
                            BC_2(loc)*B_O_Y(tr,s,3)+...
                            BC_3(loc)*B_O_Y(tr,s,5);
            Phi_temp=0;
            
            %Determine the value of Phi at location loc by summing the
            %nonzero basis function (3 vertices * 3 basis functions per
            %vertex = 9). Note that each of there basis functions has to be
            %constructed again from its (6) Bezier Ordinates. 
            %Loop over all vertices of the main triangle
            for v=Tri(tr,:)

                %loop over all basis functions of the vertex
                for q=1:3
                    %Calculating the particular basis functions from its
                    %Bezier Ordinates and adding up the contribution.
                    Phi_temp=Phi_temp+c(v,q)* (...
                        B_O(v,q,tr,s,1)*BC_1(loc)^2+...
                        B_O(v,q,tr,s,2)*2*BC_1(loc)*BC_2(loc)+...
                        B_O(v,q,tr,s,3)*BC_2(loc)^2+...
                        B_O(v,q,tr,s,4)*2*BC_2(loc)*BC_3(loc)+...
                        B_O(v,q,tr,s,5)*BC_3(loc)^2+...
                        B_O(v,q,tr,s,6)*2*BC_3(loc)*BC_1(loc) );
                end
            end
            
            %Insert Phi_temp
            Phi(tr,s,loc)=Phi_temp;
        end
%         scatter(squeeze(Phi_X(tr,s,:)),squeeze(Phi_Y(tr,s,:)))
    end
end

xlin = linspace(0,1,n_x*(resolution-1)+1);
ylin = linspace(0,1,n_y*(resolution-1)+1);

[X,Y] = meshgrid(xlin,ylin);
SI = scatteredInterpolant(reshape(Phi_X,[numel(Phi_X),1]),...
                          reshape(Phi_Y,[numel(Phi_Y),1]),...
                          reshape(Phi  ,[numel(Phi  ),1]) );

%The exact solution and the projection, and the difference between the two.
Z_exact  = f(X,Y);
Z_proj   = SI(X,Y);
Z_diff   = Z_proj-Z_exact;

figure
surf(X,Y,Z_exact,'EdgeColor','none')
title('Function f')

figure
surf(X,Y,Z_proj,'EdgeColor','none')
title('Projection of f')

figure
surf(X,Y,Z_diff,'EdgeColor','none')
title('Error')
       
clear i j k loc q s tr v xlin ylin 
toc
%% Calculate the L2 error

Error=sqrt(1/numel(X)*sum(sum(Z_diff.^2)));
fprintf('For (h_x,h_y)=(%0.5f,%0.5f), the L_2 error is %0.5e \n',...
    1/(n_x-1),1/(n_y-1),Error);

% %%
% close all
% figure
% surf(X,Y,Z_proj,'EdgeColor','none')
% xlabel('x')
% ylabel('y')
% zlabel('\phi')
% 
% figure
% trimesh(Tri,V_X,V_Y)
% hold on
% contour(X,Y,Z_proj)
% scatter(0.7,0.55,'r','filled')
% hold off
% xlabel('x')
% ylabel('y')










