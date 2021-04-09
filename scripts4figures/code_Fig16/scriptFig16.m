
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) copy the file "SCARLET_v1_GMDplotFig16.m" inside the folder
% "source_code"
% 3) set the Matlab path as 'SCARLET_v1/source_code'
% 4) Uncomment the part of the text that generates the figure you are
% interested in (uncomment between signs "...%%%%%%%%%%%%..."). A part is
% dedicated to reproduce ellipsoids; the other one to reproduce spheres.
% 5) Set "sort_boole = 0;" for sequential displacement; "sort_boole = 1;"
% for sequential one
% 6) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of plots reported in Fig. 16

close all;
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ELLIPSOIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO BE MODIFIED BY THE USER: select "0" for not sequential, "1" for
% sequential
sort_boole = 0;

% -------------------------- %
% Import the STL shapes already filled with spheres: the structure of the
% input file for this test is made by the coarse ellipsoid at the center
% acting as a core; then the coarse ellipsoid and the fine one in position
% 1 and 2 respectively.

load('../STL_examples/Fig16/aggr_27_sk21.mat');                            % Now we have agg_27_sk21 in the workspace

% Core size
core_size = 270; 

% Sizes of the coating for aggregate 27sk21 (in micrometers)
a1 = 15*ones(1,258);
a2 = 23*ones(1,138);
a3 = 32*ones(1,62);
a4 = 40*ones(1,24);
a5 = 48*ones(1,22);
a6 = 56*ones(1,8);
a7 = 65*ones(1,6);
a8 = 73*ones(1,1);
a9 = 80*ones(1,2);

% Global population made of 6 repetitions -> we discard a1 and a2 sizes for
% computational efficiency

vettore_taglie_2 = [a3 a4 a5 a6 a7 a8 a9 ...
                    a3 a4 a5 a6 a7 a8 a9 ...
                    a3 a4 a5 a6 a7 a8 a9 ...
                    a3 a4 a5 a6 a7 a8 a9 ...
                    a3 a4 a5 a6 a7 a8 a9 ...
                    a3 a4 a5 a6 a7 a8 a9];                  
    
% Parameters for ellispoid scaling (i.e. the size provided in SCARLET-1.0
% is the maximum size. The above sizes are the equivalent diameter)
alfa_fine = 1.2;
beta_fine = 2.6;
alfa_coarse = 1.14;
beta_coarse = 1.35;  



for id_ellis = 1:length(vettore_taglie_2)
       
       if vettore_taglie_2(id_ellis)<63
           vettore_taglie(id_ellis) = vettore_taglie_2(id_ellis) * (alfa_fine*beta_fine)^(1/3);
           
       else
           vettore_taglie(id_ellis) = vettore_taglie_2(id_ellis) * (alfa_coarse*beta_coarse)^(1/3);
           
       end
       
end
                
cone = 1;
rays = 1;
euler = 1; 

% TO BE MODIFIED BY THE USER: in order to just have one line set
% N_repetitions = 1 (it speeds up the code for reviewers)
N_repetitions = 5;

for i=1:N_repetitions
    [output_st, vector_collisions] = SCARLET_v1_GMDplotFig16(agg_27_sk21, vettore_taglie, ...
                                                             core_size, sort_boole, ...
                                                             cone, rays, euler);
    matrix_data(i,:) = vector_collisions;
end
  

for i=1:length(vector_collisions)
    
    y(i) = mean(matrix_data(:,i));
    dy(i) = std(matrix_data(:,i));
    
end

x = [1:1:length(vector_collisions)];

figure(1)
shadedErrorBar(x, y, dy)
ylim([0 0.9]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                               SPHERES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % TO BE MODIFIED BY THE USER: select "0" for not sequential, "1" for
% % sequential
% sort_boole = 0;
% 
% % -------------------------- %
% % Import the STL shapes for spheres (we use the same STL file of Fig.11 and
% % Fig.12
% 
% [sph_particle, fv] = fromStlToSpheres('../STL_examples/Fig11_Fig12/ParticleParticle16.stl', 1, 1, 1, 0);
% 
% cored_struct(1).fv = fv;
% cored_struct(1).sphere_struct = sph_particle;
% cored_struct(2).fv = fv;
% cored_struct(2).sphere_struct = sph_particle;
% cored_struct(3).fv = fv;
% cored_struct(3).sphere_struct = sph_particle;
% 
% % Core size
% core_size = 270; 
% 
% % Sizes of the coating for aggregate 27sk21 (in micrometers)
% a1 = 15*ones(1,258);
% a2 = 23*ones(1,138);
% a3 = 32*ones(1,62);
% a4 = 40*ones(1,24);
% a5 = 48*ones(1,22);
% a6 = 56*ones(1,8);
% a7 = 65*ones(1,6);
% a8 = 73*ones(1,1);
% a9 = 80*ones(1,2);
% 
% % Global population made of 6 repetitions -> we discard a1 and a2 sizes for
% % computational efficiency
% 
% vettore_taglie_2 = [a3 a4 a5 a6 a7 a8 a9 ...
%                     a3 a4 a5 a6 a7 a8 a9 ...
%                     a3 a4 a5 a6 a7 a8 a9 ...
%                     a3 a4 a5 a6 a7 a8 a9 ...
%                     a3 a4 a5 a6 a7 a8 a9 ...
%                     a3 a4 a5 a6 a7 a8 a9];           
%                 
% vettore_taglie = vettore_taglie_2;
% 
%                 
% cone = 1;
% rays = 1;
% euler = 1; 
% 
% % TO BE MODIFIED BY THE USER: in order to just have one line set
% % N_repetitions = 1 (it speeds up the code for reviewers)
% N_repetitions = 5;
% 
% for i=1:N_repetitions
%     [output_st, vector_collisions] = SCARLET_v1_GMDplotFig16(agg_27_sk21, vettore_taglie, ...
%                                                              core_size, sort_boole, ...
%                                                              cone, rays, euler);
       
%     matrix_data(i,:) = vector_collisions;
% end
%   
% 
% for i=1:length(vector_collisions)
%     
%     y(i) = mean(matrix_data(:,i));
%     dy(i) = std(matrix_data(:,i));
%     
% end
% 
% x = [1:1:length(vector_collisions)];
% 
% figure(1)
% shadedErrorBar(x, y, dy)
% ylim([0 0.9]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%