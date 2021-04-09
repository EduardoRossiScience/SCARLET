
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) copy the file "SCARLET_v1_GMDplotFig9new.m" inside the folder
% "source_code"
% 3) set the Matlab path as 'SCARLET_v1/source_code'
% 4) Uncomment the part of the text that generates the figure you are
% interested in (uncomment between signs "...%%%%%%%%%%%%..."
% 5) set the right number of spheres in "fromStlToSpheres"
% 6) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Fig. 9a 

close all;

% Generation of the sphere-composite representation

[sph_particle, fv] = fromStlToSpheres('../STL_examples/regular.stl', 100, 5, 0, 0);

cored_struct(1).fv = fv;
cored_struct(1).sphere_struct = sph_particle;
cored_struct(2).fv = fv;
cored_struct(2).sphere_struct = sph_particle;
N_particles = 1;
core_size = 100;
size_min = 2;
size_max = 2;
cone = 1;
rays = 1;
euler = 1;                                                                 % Always = 1 for spheres

[output_st] = SCARLET_v1_GMDplotFig9new(cored_struct, N_particles, ...
                                     core_size, size_min, size_max, ...
                                     cone, rays, euler);
output_st.packing_info.porosity  


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 9b 
% 
% close all;
% 
% % Generation of the sphere-composite representation
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/pyramid.stl', 100, 5, 0, 0);
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% N_particles = 1;
% core_size = 100;
% size_min = 2;
% size_max = 2;
% cone = 1;
% rays = 1;
% euler = 1;                                                                 % Always = 1 for spheres
% 
% [output_st] = SCARLET_v1_GMDplotFig9new(cored_struct, N_particles, ...
%                                      core_size, size_min, size_max, ...
%                                      cone, rays, euler);
% output_st.packing_info.porosity  
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 9c
% 
% close all;
% 
% % Generation of the sphere-composite representation
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/Menger2.stl', 100, 5, 0, 0);
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% N_particles = 1;
% core_size = 100;
% size_min = 2;
% size_max = 2;
% cone = 1;
% rays = 1;
% euler = 1;                                                                 % Always = 1 for spheres
% 
% [output_st] = SCARLET_v1_GMDplotFig9new(cored_struct, N_particles, ...
%                                      core_size, size_min, size_max, ...
%                                      cone, rays, euler);
% output_st.packing_info.porosity  
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 9d
% 
% close all;
% 
% % Generation of the sphere-composite representation
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/Menger3.stl', 100, 5, 0, 0);
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% N_particles = 1;
% core_size = 100;
% size_min = 2;
% size_max = 2;
% cone = 1;
% rays = 1;
% euler = 1;                                                                 % Always = 1 for spheres
% 
% [output_st] = SCARLET_v1_GMDplotFig9new(cored_struct, N_particles, ...
%                                      core_size, size_min, size_max, ...
%                                      cone, rays, euler);
% output_st.packing_info.porosity  
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 9e
% 
% close all;
% 
% % Generation of the sphere-composite representation
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/L_shaped.stl', 100, 5, 0, 0);
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% N_particles = 1;
% core_size = 100;
% size_min = 2;
% size_max = 2;
% cone = 1;
% rays = 1;
% euler = 1;                                                                 % Always = 1 for spheres
% 
% [output_st] = SCARLET_v1_GMDplotFig9new(cored_struct, N_particles, ...
%                                      core_size, size_min, size_max, ...
%                                      cone, rays, euler);
% output_st.packing_info.porosity 

