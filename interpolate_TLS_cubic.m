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

function [value_gp] = interpolate_TLS_cubic(particles_X,particles_Y,...
    value_p,volume_p,GaussPoints_X,GaussPoints_Y,w_gp,...
    triangles_active_boundary,triangles_inactive, triangles,...
    conservation)

%The number of triangles, particles and Gauss points
n_triangles=size(triangles.ConnectivityList,1);
n_gp = length(GaussPoints_X);

value_gp=zeros(n_gp,1);

%The triangles and barycentric coordinates of all particles and positions
%of interest
[p_tr,BC_P] = pointLocation(triangles,particles_X,particles_Y);
[gp_tr,BC_gp] = pointLocation(triangles,GaussPoints_X,GaussPoints_Y);

%Loop over the triangles
for tri = 1:n_triangles
    if (ismember(tri,triangles_inactive) || ...
            ismemeber(tri,triangles_active_boundary))
        continue    %skip tri if tri is an inactive/boundary triangle
    end
    %Particles and Guass points in the considered triangle
    p_tri = find(p_tr==tri);
    gp_tri = find(gp_tr==tri);
    
    %Barycentric coordinates of particles
    eta1=BC_P(p_tri,1)';
    eta2=BC_P(p_tri,2)';
%     eta3=BC_P(p_tri,3)';
    
    %Barycentric coordinates of Gauss Points
    eta1_gp=BC_gp(gp_tri,1)';
    eta2_gp=BC_gp(gp_tri,2)';
%     eta3_gp=BC_gp(gp_tri,3)';
    
    %The integral of BC eta_i^n over the triangle is equal to  
    %1/(1/2)*\int_\Omega \eta_i^n d\Omega = 
    %2*\int_0^1 (1-\eta_i)*\eta_i^n d\eta_i =
    %2*\int_0^1 \eta_i^n-\eta_i^(n+1) d\eta_i =
    %2*\left[1/(n+1)*\eta_i^(n+1)-1/(n+2)*\eta_i^(n+2)\right]_0^1 =
    %2*(1/(n+1)-1/(n+2)) = 2/((n+1)*(n+2))
        
    %The basis functions
    phi = [ones(size(eta1));    eta1-1/3;           eta2-1/3; ...
           eta1.^2-1/6;         eta1.*eta2-1/12;    eta2.^2-1/6;    ...
           eta1.^3-1/10;        eta1.^2.*eta2-1/30; ...
           eta1.*eta2.^2-1/30;  eta2.^3-1/10]';
    
    %A and b in Ac=b, for determining the least squares approximation
    A = phi'*(volume_p(p_tri)*ones(1,size(phi,2)).*phi);
    b = [   ones(size(eta1));    eta1 - 1/3;        eta2 - 1/3; ...
            eta1.^2-1/6;         eta1.*eta2-1/12;   eta2.^2-1/6;    ...
            eta1.^3-1/10;        eta1.^2.*eta2-1/30; ...
            eta1.*eta2.^2-1/30;  eta2.^3-1/10      ]*...
        (volume_p(p_tri).*value_p(p_tri));
    
    %Impose that the first coefficient is equal to the average over the
    %triangle for mass conservation
    if conservation == 1
        A(1,:)=0;
        b(1)=volume_p(p_tri)'*value_p(p_tri)/sum(w_gp(gp_tri));
    %     b(1)=sum(value_p(p_tri))/numel(p_tri);
        b(2:end)=b(2:end)-b(1)*A(2:end,1);
        A(:,1)=0;
        A(1,1)=1;
    end
    
    %Find the coefficients for the reconstruction
    c = A\b;
    
    %Find the interpolated values in the Gauss points
    phi_gp = [ones(size(eta1_gp)); eta1_gp-1/3;         eta2_gp-1/3;...
              eta1_gp.^2-1/6;  eta1_gp.*eta2_gp-1/12;   eta2_gp.^2-1/6;...
              eta1_gp.^3-1/10; eta1_gp.^2.*eta2_gp-1/30; ...
              eta1_gp.*eta2_gp.^2-1/30;  eta2_gp.^3-1/10];
    value_gp(gp_tri) = phi_gp'*c;
end

end
    

