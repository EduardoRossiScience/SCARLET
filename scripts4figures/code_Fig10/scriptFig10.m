
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) copy the file "SCARLET_v1_GMDplotFig10.m" inside the folder
% "source_code"
% 3) set the Matlab path as 'SCARLET_v1/source_code'
% 4) Uncomment the part of the text that generates the figure you are
% interested in (uncomment between signs "...%%%%%%%%%%%%..."
% 5) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of points associated with spheres in Fig. 10 (black dots in the plot) 

close all;

% Generation of the sphere

[sph_particle, fv] = fromStlToSpheres('../STL_examples/sphere.stl', 1, 10, 1, 0);

cored_struct(1).fv = fv;
cored_struct(1).sphere_struct = sph_particle;
cored_struct(2).fv = fv;
cored_struct(2).sphere_struct = sph_particle;
N_particles = 100;
core_size = 10;
size_min = 10;
size_max = 10;
cone = 90;
rays = 10;
euler = 1;                                                                 % Always = 1 for spheres

[output_st] = SCARLET_v1_GMDplotFig10(cored_struct, N_particles, ...
                                     core_size, size_min, size_max, ...
                                     cone, rays, euler);

x_coordinate_plot = core_size/output_st.packing_info.size                 % x coordinate of the point in Fig.9 (spheres)
y_coordinate_plot = 1-output_st.packing_info.porosity                     % y coordinate of the point in Fig.9 (spheres)

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of points associated with ellipsoids in Fig. 9 (black dots in the plot) 
% 
% close all;
% 
% % Generation of the ellipsoid
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/ellipse.stl', 300, 1, 1, 0);
% 
% 
% % IMPORTANT: modify the following three numbers according to the ones
% % reported in the round brackets of Fig.10 
% cone = 10;                                                                 
% rays = 1;
% euler = 1;
% 
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% N_particles = 50;                                                          % For points with gamma<0.1 set "N_particles = 100";
% core_size = 100;
% size_min = 100;
% size_max = 100;
% 
% [output_st] = SCARLET_v1_GMDplotFig10(cored_struct, N_particles, ...
%                                      core_size, size_min, size_max, ...
%                                      cone, rays, euler);
% 
% x_coordinate_plot = core_size/output_st.packing_info.size                 % x coordinate of the point in Fig.9 (ellipsoids)
% y_coordinate_plot = 1-output_st.packing_info.porosity                     % y coordinate of the point in Fig.9 (ellipsoids)
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%