
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) copy the file "SCARLET_v1_GMDplotFig12_13.m" inside the folder
% "source_code"
% 3) set the Matlab path as 'SCARLET_v1/source_code'
% 4) Uncomment the part of the text that generates the figure you are
% interested in (uncomment between signs "...%%%%%%%%%%%%..."
% 5) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of points associated with Fig. 12 and Fig.13. 
% Please note that the present script just generates one of the 16 plots
% reported in each of the two figures (Fig.12 and Fig.13). In order to have
% all the 32 figures the user must change the ellipsoids and the setup for
% the simulation according to the ones reported in the text

close all;
clear all; 

% -------------------------- %
% TO BE MODIFIED BY THE USER 
% Generation of the ellispoid: here the external user must change the final
% name of the STL file used (i.e. from "ParticleParticle1.stl" to
% "ParticleParticle16.stl" depending on the figure.
% Warning! Do not modify the path (i.e. ellipsoids are inside the folder
% "Fig12_Fig13")

[sph_particle, fv] = fromStlToSpheres('../STL_examples/Fig12_Fig13/ParticleParticle16.stl', 1, 1, 1, 0);

% TO BE MODIFIED BY THE USER
% This is the setup that makes the difference between Fig.12 and Fig.13.
% Computations made for Fig.12 are usually more demanding and they
% generally require much more time (and a cluster9
cone = 1;
rays = 1;
euler = 1; 

% -------------------------- %

cored_struct(1).fv = fv;
cored_struct(1).sphere_struct = sph_particle;
cored_struct(2).fv = fv;
cored_struct(2).sphere_struct = sph_particle;

ratios = [1 2 3 4 5 6 7 8 9 10];

% coating size
coating_size = 10;

N_particles = 1;
N_repetitions = 5;
porosity = zeros(N_repetitions, length(ratios));

for i = 1:N_repetitions
    
    for j = 1:length(ratios)
    
        % core size
        core_size = coating_size*ratios(j);
    
        [output_st] = SCARLET_v1_GMDplotFig12_13(cored_struct, N_particles, ...
                                              core_size, coating_size, coating_size, ...
                                              cone, rays, euler);
                                          
        porosity(i,j) = output_st.packing_info.porosity;
    
    end
end

for i=1:length(ratios)
    
    avg_plot(i) = mean(porosity(:,i));
    std_plot(i) = std(porosity(:,i));
    
end

% figure(1)
% errorbar(ratios, avg_plot, std_plot, '.')

figure(1)
shadedErrorBar(ratios, avg_plot, std_plot)


