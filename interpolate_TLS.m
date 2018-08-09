%Interpolate a function known in particle points (particles_X,particles_Y)
%with values value_p to get the values value_gp in the Gausspoints 
%(GaussPoints_X,GaussPoints_Y). The interpolation uses a discontinuous 
%Taylor Least Squares reconstruction over each element in triangles. Over
%each triangle, the function will be as a 2D parabola. Over the whole
%domain, discontinuities will appear over element edges. For integration
%purposes, this is not a problem. The TLS reconstruction is described for
%1D by E. Wobbes, and more about this can be found in a report of hers,
%"Conservative Taylor Least Squares reconstruction for material point 
%methods"
%In this script, particles_X, particles_Y and value_p are all vectors of
%length n_particles. GaussPoints_X, GaussPoints_Y and value_gp will be
%vectors of length n_GaussParticles and triangles is a triangulation
%containing all the particles and Gauss points. 

function [value_gp] = interpolate_TLS(particles_X,particles_Y,value_p,...
    volume_p,GaussPoints_X,GaussPoints_Y,w_gp,triangles_active_boundary, ...
                triangles_inactive, triangles,mass_conservation)

%The number of triangles, particles and Gauss points
n_triangles=size(triangles.ConnectivityList,1);
% n_particles=length(particles_X);
n_gp = length(GaussPoints_X);

value_gp=zeros(n_gp,1);

% %The vertices of each triangle
% V1=triangles.Points(triangles.ConnectivityList(:,1),:);
% V2=triangles.Points(triangles.ConnectivityList(:,2),:);
% V3=triangles.Points(triangles.ConnectivityList(:,3),:);
% 
% %The surface of each triangle
% sur = abs(V1(:,1).*(V2(:,2)-V3(:,2))+...
%          V2(:,1).*(V3(:,2)-V1(:,2))+...
%          V3(:,1).*(V1(:,2)-V2(:,2)) )/2;

%The triangles and barycentric coordinates of all particles and positions
%of interest
[p_tr,BC_P] = pointLocation(triangles,particles_X,particles_Y);
[gp_tr,BC_gp] = pointLocation(triangles,GaussPoints_X,GaussPoints_Y);

%Loop over the triangles
for tri = 1:n_triangles
    if (ismember(tri,triangles_inactive) || ...
            ismember(tri,triangles_active_boundary))
        continue    %skip tri if tri is an inactive/boundary triangle
    end
    %Particles and Guass points in the considered triangle
    p_tri = find(p_tr==tri);
    gp_tri = find(gp_tr==tri);
    
    %Barycentric coordinates of particles
    eta1=BC_P(p_tri,1)';
    eta2=BC_P(p_tri,2)';
    eta3=BC_P(p_tri,3)';
    
    %Barycentric coordinates of Gauss Points
    eta1_gp=BC_gp(gp_tri,1)';
    eta2_gp=BC_gp(gp_tri,2)';
    eta3_gp=BC_gp(gp_tri,3)';
    
    %The basis functions
    phi = [ones(size(eta1)); eta1-1/3; eta2-1/3; eta1.^2-1/6;...
            eta2.^2-1/6; eta3.^2-1/6]';
    
    %A and b in Ac=b, for determining the least squares approximation
%     A = phi'*phi;
%     b = [ones(size(eta1)); eta1 - 1/3; eta2 - 1/3; eta1.^2-1/6; ...
%         eta2.^2-1/6; eta3.^2-1/6] * ...
%         value_p(p_tri);
    A = phi'*(volume_p(p_tri)*ones(1,size(phi,2)).*phi);
    b = [ones(size(eta1)); eta1 - 1/3; eta2 - 1/3; eta1.^2-1/6; ...
        eta2.^2-1/6; eta3.^2-1/6] * ...
        (volume_p(p_tri).*value_p(p_tri));
    
    if (mass_conservation==1)
        %Impose that the first coefficient is equal to the average over the
        %triangle for mass conservation
        A(1,:)=0;
        b(1)=volume_p(p_tri)'*value_p(p_tri)/sum(w_gp(gp_tri));
%         b(1)=sum(value_p(p_tri))/numel(p_tri);
        b(2:end)=b(2:end)-b(1)*A(2:end,1);
        A(:,1)=0;
        A(1,1)=1;
    end 
    
    if rcond(A)<1e-16;
    end
    
    %Find the coefficients for the reconstruction
    c = A\b;
    
    %Find the interpolated values in the Gauss points
    phi_gp = [ones(size(eta1_gp)); eta1_gp-1/3; eta2_gp-1/3;...
              eta1_gp.^2-1/6; eta2_gp.^2-1/6; eta3_gp.^2-1/6];
    value_gp(gp_tri) = phi_gp'*c;
end

end
    

