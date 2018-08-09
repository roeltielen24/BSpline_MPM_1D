%Calculate the Cartesian coordinates of the Gauss points for each triangle
%in the given triangulation. The result for the x-locations is a matrix of
%size (n_triangles, n_GaussPointsPerTriangle), i.e. the locations for the 
%Gauss point gp in triangle i is found in GaussPoints_X/Y(i,gp), and each
%point gp in triangle i has weight Gauss_weight(i,gp)
function [GaussPoints_X,GaussPoints_Y,Gauss_weight]=...
            Gauss_Points(Tri)

Gauss_degree=5;
        
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

%Order the points of the triangles such that in the rows are the triangles
%and in the columns the vertices 1 to 3 of the that triangle.
V_X=reshape(Tri.Points(Tri.ConnectivityList,1),[],3);
V_Y=reshape(Tri.Points(Tri.ConnectivityList,2),[],3);

GaussPoints_X=V_X*[Gauss_t1;Gauss_t2;Gauss_t3];
GaussPoints_Y=V_Y*[Gauss_t1;Gauss_t2;Gauss_t3];
Tri_area=1/2*abs((V_X(:,2)-V_X(:,1)).*(V_Y(:,3)-V_Y(:,1))-...
                 (V_X(:,3)-V_X(:,1)).*(V_Y(:,2)-V_Y(:,1)) );
Gauss_weight=Tri_area*Gauss_W;
end