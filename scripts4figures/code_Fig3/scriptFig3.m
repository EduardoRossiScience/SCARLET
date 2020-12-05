
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) set the Matlab path as 'SCARLET_v1/source_code'
% 3) Uncomment the part of the text that generates the figureyou are
% interested in
% 4) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Fig. 3a, 3b
close all;

[spheres_struct, fv] = fromStlToSpheres('../STL_examples/particle.stl', 300, 10, 1, 0);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 3c
% close all;
% 
% [spheres_struct, fv] = fromStlToSpheres('../STL_examples/particle.stl', 300, 10, 1, 1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 3d, 3e
% close all;
% 
% [spheres_struct, fv] = fromStlToSpheres('../STL_examples/particle2.stl', 300, 10, 1, 0);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generation of Fig. 3f
% close all;
% 
% [spheres_struct, fv] = fromStlToSpheres('../STL_examples/particle2.stl', 300, 10, 1, 1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%