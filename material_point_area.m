%From the set of material points (a.k.a. particles), determine the areas
%corresponding to each of those particle points by using a voronoi
%decomposition (i.e. each point in the domain is in the area of the
%particle it is closest to). These areas are stored by using two large
%vectors for the x-positions and the y-positions of all vertices of these
%areas and an index vector, determining which of these vertices are of the
%area of each material point. Finally, the surface of each of these areas
%is also returned, which will be used as the weight of each particle. 

function [area_vertex_X,area_vertex_Y,area_vertex_indices,area_weight] =...
    material_point_area(constant,particles_X,particles_Y)

%Mirror the particles in the boundaries to form the boundary for each
%surface
particles_temp_X=[particles_X;2*constant.length-particles_X;...
                  particles_X; -particles_X; particles_X];
particles_temp_Y=[particles_Y; particles_Y; 2*constant.width-particles_Y;
                  particles_Y; -particles_Y];

% %plot the material points, the mirrored points and the corresponding
% %surfaces
% figure
% voronoi(particles_temp_X,particles_temp_Y)
% title('Material points and weigth surfaces')

%Construct the polygons of all the surfaces of the particle points, where
%V_temp are the vertices of the polygons and voronoi_cells_temp the cells
%with the information which vertices make up the polygon.
[area_vertex_temp,area_vertex_indices_temp]=...
    voronoin([particles_temp_X,particles_temp_Y]);

%Find the list of polygon vertices that are outside the original domain.
list_out=find( (area_vertex_temp(:,1)<-1e-15) +...
               (area_vertex_temp(:,1)>constant.length+1e-15) +...
               (area_vertex_temp(:,2)<-1e-15) + ...
               (area_vertex_temp(:,2)>constant.width+1e-15) );
list_in=(1:size(area_vertex_temp,1))';
list_in(list_out)=[];          

area_vertex_X=area_vertex_temp(list_in,1);
area_vertex_Y=area_vertex_temp(list_in,2);

%All the vertices outside the original domain will not be used anymore, and
%are therefore discarded. New indices shall therefore be given. index
%change contains the new index of each of the old indices. The vertices
%outside the domain shall receive index 0, the rest will get index 1,2,3...
index_change=1:size(area_vertex_temp,1);
index_change(list_out)=0;
index_change(list_in)=1:length(list_in);

area_vertex_indices=cell(1,length(particles_X));

area_weight=ones(length(particles_X),1);

%Filter the polygons of the original material points from the list of all
%polygons
num_active_cells=0;
for i=1:numel(area_vertex_indices_temp)
    cell_as_matrix=cell2mat(area_vertex_indices_temp(i));
    %if none of the vertices of a polygon are in the list with outlying
    %vertices, then the full polygon is in the original domain.
    if sum(ismember(cell_as_matrix,list_out))==0
        %shift the indices to include only the relevant vertices
        cell_as_matrix=index_change(cell_as_matrix);
        num_active_cells=num_active_cells+1;
        area_vertex_indices(num_active_cells)={cell_as_matrix};
        [~,area_weight(num_active_cells)]=convhull(...
            area_vertex_X(cell_as_matrix),area_vertex_Y(cell_as_matrix) );
    end

end
area_vertex_indices=area_vertex_indices(1:num_active_cells);

end



