function [x sol vel] = exact_solution(rho, E, p0, g, H, x0,t )
% x0 = 1; % original node position
% t=0:0.01:0.5;
% 
% 
% rho = 10^3; %density
% E = 10^5; %stiffness
% p0 = -50*10^3; %load
% g = 9.8;
% H = 1; %height 

%% Analytical solution for displacement
%calculate coefficients 
for n=1:20
    u(n) = (8*H*(2*pi*p0*n*(-1)^n + 2*rho*g*H - pi*p0*(-1)^n))/...
    ((4*n^2 - 4*n +1)*(2*n-1)*(pi^3)*E);
end

%Calculation of the analytical solution for displacement
C_x0 = (rho*g*x0^2)/(2*E) + (p0-rho*g*H)*x0/E;
for l=1:length(t)
    for k=1:20
        part(k)=u(k)*cos(sqrt(E/rho)*(2*k-1)*pi*t(l)/(2*H))*sin((2*k-1)*...
            pi*x0/(2*H));
        vel_part(k) = -u(k)*sqrt(E/rho)*(2*k-1)*pi/(2*H)*... 
            sin(sqrt(E/rho)*(2*k-1)*pi*t(l)/(2*H))*sin((2*k-1)*pi*x0/(2*H));  

    end
    sol(l) = C_x0 + sum(part);
    vel(l) = sum(vel_part); 
    clear part vel_part;
end

%Find the position
x = x0+sol;
end