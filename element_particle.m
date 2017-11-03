function [n] = element_particle(A)
%provided the row of the elements_particles matrix corresponding with the 
%index of an element, it returns the indices of the particles
%contained in that element
[m,n] = find(A);