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
close all

%% Parameters (This part, and only this part, is for tweaking!)
%Computation time measurement
tic
disp('Initialising...')

%Dimension of the problem, do not change this, the code is specialised for
%2D!
dim=2;

%Enable lumping
lumping = 1;

%Function to be approximated
f=@(x,y) 0.1*sin(pi.*x);
% f=@(x,y) 0.1*sin(pi.*x)*sin(pi.*y);
% f=@(x,y) 1.*x;
% f=@(x,y) heaviside(x-0.7).*(x-0.7).*heaviside(y-0.7).*(y-0.7);

%In case we want to impose a essential homogeneous boundary condition on
%all sides in correspondence with the approximated function, put BC
%(Boundary Condition) to 1.
BC=2;

%Number of grid points in each direction (at least 2). 
n_x=16+1;
n_y=16+1;
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
%% Define the piece-wise linear Lagrange functions
tic
disp('Calculating Bezier Ordinates of the basis functions...')

toc
%% Construct the Mass matrix M and the Right Hand Side (RHS)
tic
disp('Constructing Mass matrix and Right Hand Side...')
%The figures and equations refered to can be found in  
%"On calculating normalized Powell-Sabin B-splines", by Paul Dierckx.

%Preallocate the mass matrix. M_ij contains \int phi_i*phi_j dOmega, where
%i=3*(v-1)+q, with v the vertex number (from 1 to Num_Ver) and q the number
%of the vertex of the PS triangle (from 1 to 3).
M=zeros(Num_Ver);

%Preallocate the RHS
RHS=zeros(Num_Ver,1);

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

Gauss_BC=[Gauss_t1;...
          Gauss_t2;...
          Gauss_t3];

%loop over all i (rows) of M_ij and RHS_i
%loop over all triangles
for tri = 1:Num_Tri
    ver_tri = Tri.ConnectivityList(tri,:);
    V=[V_X(ver_tri)';V_Y(ver_tri)'];
    area_tri = abs( V(1,1)*(V(2,2)-V(2,3)) + ...
                    V(1,2)*(V(2,3)-V(2,1)) + ...
                    V(1,3)*(V(2,1)-V(2,2)) )/2;
    for i=1:3
        for j=1:3
            M(ver_tri(i),ver_tri(j)) = M(ver_tri(i),ver_tri(j)) + ...
                area_tri * Gauss_BC(i,:)*Gauss_BC(j,:)';
        end
        x_gp = V(1,:)*Gauss_BC;
        y_gp = V(2,:)*Gauss_BC;
        RHS(ver_tri(i)) = RHS(ver_tri(i)) + ...
            area_tri * Gauss_BC(i,:)*f(x_gp',y_gp');
    end
    
end

clear gp gp_X gp_Y i j M_temp qi qj s tr Tri_vi Tri_vivj Tri_vj ...
    v1 v2 v3 vi vj    
toc


%% Implement boundary conditions

%In case of a homogeneous Dirichlet boundary condition on all sides
if BC==1
    Bound_south=find(V_Y == 0);
    Bound_north=find(V_Y == 1);
    Bound_west= find(V_X == 0);
    Bound_east= find(V_X == 1);
    
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
    Bound_west= find(V_X == 0);
    Bound_east= find(V_X == 1);
    
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
else
    c=M\RHS;
end

% c=zeros(Num_Ver,1); c(1)=1;

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
Phi_X=zeros(Num_Tri,resolution*(resolution+1)/2);
Phi_Y=zeros(Num_Tri,resolution*(resolution+1)/2);
Phi  =zeros(Num_Tri,resolution*(resolution+1)/2);

%Loop over all main triangles
for tr=1:Num_Tri     
    %loop over all locations in the subtriangle
    for loc=1:length(BC_1)
        ver_tri = Tri.ConnectivityList(tr,:);
        %Determine the X and Y values of the current location
        Phi_X(tr,loc)=BC_1(loc)*V_X(ver_tri(1))+...
                      BC_2(loc)*V_X(ver_tri(2))+...
                      BC_3(loc)*V_X(ver_tri(3));
        Phi_Y(tr,loc)=BC_1(loc)*V_Y(ver_tri(1))+...
                      BC_2(loc)*V_Y(ver_tri(2))+...
                      BC_3(loc)*V_Y(ver_tri(3));
        Phi(tr,loc) = c(ver_tri)'*[BC_1(loc);BC_2(loc);BC_3(loc)];
    end
%     scatter(squeeze(Phi_X(tr,s,:)),squeeze(Phi_Y(tr,s,:)))
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












