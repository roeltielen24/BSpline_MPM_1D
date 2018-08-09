close all; clear; %clc;

constant = struct('dim',2,'density',1,'E',100,'g',0,'load',0,...
                  'length',pi/2,'width',pi/2,'alpha',0,'v_0',0.01,...
                  'Poisson_ratio',0);

% f = @(x,y) x.^2.*y;
f = @(x,y) 1./(x+1).^2.*sin(y);
% f = @(x,y) (x_p+y_p>1).*sin(x_p).*sin(y_p);

%m is the number of particles in each direction
m=10;
x_p=reshape(linspace(0.04,constant.length-0.06,m+1)'*ones(1,m+1),[],1);
y_p=reshape(ones(m+1,1) *linspace(0.05,constant.width-0.05,m+1),[],1);
z_p = f(x_p,y_p);

[~,~,~,volume_p] = material_point_area(constant,x_p,y_p);

% n=50;
% x_gp=constant.length*rand(n,1);
% y_gp=constant.width *rand(n,1);
% x_gp=reshape(linspace(0.04,constant.length-0.06,m+1)'*ones(1,m+1),[],1);
% y_gp=reshape(ones(m+1,1) *linspace(0.05,constant.width-0.05,m+1),[],1);

Tri = delaunayTriangulation(constant.length*[0 0 1 1]',...
                            constant.width *[0 1 0 1]');

[x_gp,y_gp,w_gp] = Gauss_Points(Tri);

x_gp=reshape(x_gp',[],1);
y_gp=reshape(y_gp',[],1);
w_gp=reshape(w_gp',[],1);

figure
hold on
triplot(Tri)
scatter(x_p,y_p,'b')
scatter(x_gp,y_gp,'r')

conserve = 1;
z_gp = interpolate_TLS(x_p,y_p,z_p,volume_p,x_gp,y_gp,w_gp,Tri,conserve);
% z_gp = interpolate_TLS_Version2(x_p,y_p,z_p,volume_p,x_gp,y_gp,w_gp,Tri);
% z_gp = interpolate_TLS_Version3(x_p,y_p,z_p,volume_p,x_gp,y_gp,w_gp,Tri);

[X,Y]=meshgrid(0:constant.length/100:constant.length,...
               0:constant.width /100:constant.width);
Z = f(X,Y);

figure
hold on
surf(X,Y,Z,'linestyle','none')
scatter3(x_p,y_p,z_p,'b')
scatter3(x_gp,y_gp,z_gp,'r')

L2_error = sqrt( w_gp'*(z_gp-f(x_gp,y_gp)).^2 );

fprintf('L2-error = %f \n',L2_error);
