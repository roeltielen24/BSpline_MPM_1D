function [m] = particle_element(A)
%provided the column of the elements_particles matrix corresponding with 
%the index of a particle, it returns the index of the element
%containing that particle
[m,n] = find(A);