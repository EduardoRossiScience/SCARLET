function [output_st, vettore_poro] = SCARLET_v1_GMDplotFig17(cored_struct, vettore_taglie, ...
                                                             core_size, sort_boole, ...
                                                             cone, rays, euler)

    % INFO about the input parameters: 
    
    % "cored_struct" is the initial structure that contains all the
    % information about the shapes that will be used by the program to
    % create the aggregates. 
    %
    % It has the following inner structure:
    % - cored_struct(1).fv: faces and vertices of the polygons used for the
    % stl file.
    % - cored_struct(1).sphere_struct: the coordinates, the radius, the
    % volume and the unique id of each of the spheres used to describe the
    % closed surface contained in the ".fv" file above.
    %
    % "cored_struct" has at least two objects, i.e. the following structure:
    % cored_struct(1).fv;
    % cored_struct(1).cored_struct;
    % cored_struct(2).fv;
    % cored_struct(2).cored_struct;
    % The first object is used as the core of the aggregate. You can add as
    % many objects as you wish.
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                 DEFINE INITIAL CONDITIONS                   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Core size (in micrometers)
    closet.core_size_mu = core_size;                                             % Core size
    
    % Density of the core (kg/m^3)
    closet.core_density = 2500;  
    
    % Coating particles
    closet.N_particles = length(vettore_taglie);                           % Number of the coating particles (core not included)
                
    % Sequential displacement or not
    if (sort_boole)<1             
        closet.vector_coat_sizes_mu = vettore_taglie(randperm(length(vettore_taglie)));
    else
        closet.vector_coat_sizes_mu = sort(vettore_taglie,'descend');
    end
    
    closet.displace_in_order = 0;                                          % Boolean variable to decide if the coating is random (0) or not (1)
    
    % Cone of inspection features
    closet.cone_aperture_degree = cone;                                    % Cone aperture of inspection: [deg]
    closet.N_raysXSA = rays;                                               % Number of rays per each solid angle
    
    % Euler angles 
    closet.N_Euler_triplets = euler;                                       % Number of angles per fixed direction within the solid angle
    
    % Versus of placing particles
    closet.origin_in_the_CM = 1;                                           % Center of the cone placed on the center of mass of the inner core or not
    
    % Plotting features
    closet.alpha_transparency = 1;                                         % Transparency of figures
    closet.view_angle_1 = 35;                                              % View angle #1
    closet.view_angle_2 = 42;                                              % View angle #2
    
    % Controls over the iterations
    closet.translation_iter_max = 1000;                                    % Maximum number of iterations
    closet.translation_vec_modulus = 1;                                    % Spatial step
    closet.delta = 0.01;                                                  % Constant for the inward movement (rough)
    closet.delta_2 = 0.0005;                                               % Constant for the outward movement (fine)
    
    % Variables to speed-up the code
    closet.n_points4Feret = 1000;                                          % Number of points used to calculate the Feret diameter: 
                                                                           % all objects MUST contain at least this number of points
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                    END INITIAL CONDITIONS                   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         
                            %%% MAIN: begins %%%
                            
    warning off;
    close all;
    welcomeMessage;                        
                            
    % Main function of the program: "closet" is the structure where all the
    % input parameters and auxiliary variables are contained.
    [front_points, g_Object, closet, output] = mainFunction(closet, cored_struct);
                    
    output_st.packing_info = output;
    output_st.aggregate = g_Object;
    vettore_poro = closet.vettore_porosita;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% External functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [front_points, g_Object, closet, output] = mainFunction(closet, cored_struct)

    % Additionl routine to complete the initial conditions and to create the 
    % important structure "plusier_STLs", which will be used to not
    % overwrite the input structure "cored_struct"
    [plusier_STLs, closet] = elaborateInitialStructure(cored_struct, closet);
    
    % Initialization of empty matrices and vectors used in the code
    [closet]  = initializeEmptyMatrices(closet);
    
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%% Locating objects from the outside %%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        
     % Scaling of the core to the dimensions given by the user 
     plusier_STLs.core = realObjectCreator(plusier_STLs.core, closet.core_size_mu);
        
     % We dedicate a specific structure to the central object: c_Object
     c_Object = createCObject(plusier_STLs.core.fv, plusier_STLs.core.sphere_struct);
        
     % We create a new structure that will contain all the information
     % of all the spheres contained in all the stl files. It is required
     % for the intersection
     s_Blob = initializeBlob(c_Object);
        
     % The structure "g_Object" stands for "global Object" and it is
     % where the final aggregate is stored
     g_Object(1).s_Object = c_Object;
        
     tic
        
     %%%%%%%%%%%%%%%%%%%%%%%%%      LOOP 1      %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % Loop over each single particle to be placed around the inner core   closet.N_particles
     for id_particle = 1:closet.N_particles
        
         fprintf('Now adding particle number %i of %i to the central core\n', id_particle, closet.N_particles);
        
         % Random direction of the central axis of the cone
         versor = generateUnitVector(1);
         cone_central_dir = [versor(1); versor(2); versor(3)];
 
        
         % Table (with dimensions [closet.N_raysXSA x 3]) with all the rays 
         % within the solid angle for a single particle
           
         rays_vec = solidAngleRndPoints(closet.cone_aperture_degree, ...
                                        cone_central_dir, ...
                                        closet.N_raysXSA);
       
         % Inizialisation of the distance vector between the point on the
         % surface of the blob and the center of mass of the closest object
         % respect to the new one that must be placed
         d_CM = zeros(1,id_particle);
        
            
         %%%%%%%%%%%%%%%%%%%%%%%      LOOP 2     %%%%%%%%%%%%%%%%%%%%%%%%%%
            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
         %       LOOP OVER ALL THE RAYS CONTAINED IN THE SOLID ANGLE      %
         %                 (for one single coating particle)              %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
         % Loop over all the different directions within the solid angle:
         % in this loop we place all the different objects within the
         % selected cone (this loop is not parallelized)

         for id_single_ray = 1:closet.N_raysXSA
            
             fprintf('Now investigating ray N.%i of %i\n', id_single_ray, closet.N_raysXSA);
             
             x_dir = rays_vec(1, id_single_ray);                           % x versor for the single ray
             y_dir = rays_vec(2, id_single_ray);                           % y versor for the single ray
             z_dir = rays_vec(3, id_single_ray);                           % z versor for the single ray
                  
             % Choose if the origin of the cone is set to the center of mass 
             % of the core or on a coating particle
             
             if closet.origin_in_the_CM >0
                % origin set in the center of mass of the core
                
                origin_of_the_ray(1) = c_Object.CM(1);                     % x coordinate of the origin of the ray
                origin_of_the_ray(2) = c_Object.CM(2);                     % y coordinate of the origin of the ray
                origin_of_the_ray(3) = c_Object.CM(3);                     % z coordinate of the origin of the ray
             
             else
                % origin set on a coating particle 
                 
                vec_prov = randperm(s_Blob.n_tot_spheres);
                id_start_sphere = vec_prov(1);
                origin_of_the_ray(1) = s_Blob.sphere_struct.x_c(id_start_sphere);            % x coordinate of the origin of the ray
                origin_of_the_ray(2) = s_Blob.sphere_struct.y_c(id_start_sphere);            % y coordinate of the origin of the ray
                origin_of_the_ray(3) = s_Blob.sphere_struct.z_c(id_start_sphere);            % z coordinate of the origin of the ray
             
             end
                
             % Definition of the starting point over the BLOB surface
             % where to start the object displacement algorithm (along
             % one single ray)
             
             [front_point_C_Obj, ~] = intersectionLineSphere(s_Blob, ...
                                                             origin_of_the_ray(1), x_dir, ...
                                                             origin_of_the_ray(2), y_dir, ...
                                                             origin_of_the_ray(3), z_dir);
            
             % If the center of mass is not located inside one of the
             % spheres "front_point_C_Obj.x" can be empty
             
             if isempty(front_point_C_Obj.x)                                          
                                            
                id = IDdistCalculatorVector(origin_of_the_ray(1), origin_of_the_ray(2), origin_of_the_ray(3), ...
                                            s_Blob.sphere_struct.x_c, s_Blob.sphere_struct.y_c, s_Blob.sphere_struct.z_c);                        
                    
                origin_of_the_ray(1) = s_Blob.sphere_struct.x_c(id);
                origin_of_the_ray(2) = s_Blob.sphere_struct.y_c(id);
                origin_of_the_ray(3) = s_Blob.sphere_struct.z_c(id);                       
                                            
                [front_point_C_Obj, ~] = intersectionLineSphere(s_Blob, ...
                                                                origin_of_the_ray(1), x_dir, ...
                                                                origin_of_the_ray(2), y_dir, ...
                                                                origin_of_the_ray(3), z_dir);
             end                                           
                                                        
             % Starting point over the surface of the blob object
             start_p.x = front_point_C_Obj.x;
             start_p.y = front_point_C_Obj.y;
             start_p.z = front_point_C_Obj.z;
            
                
             %%%%%%%%%%%%%%%%%%%      LOOP 3     %%%%%%%%%%%%%%%%%%%%%%%
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
             %%%%%            LOOP OVER ALL THE EULER ANGLES          %%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             % Loop over all the possible Euler rotations of the object for
             % a given - fixed - ray.
             % This is the unique parallel loop of the code: the
             % parallelization is done over each triplet of the Euler's
             % angles.
            
                
             % -------------- Scaling of the coating: BEGIN ------------
             % Here we choose randomly a shape from the list of objects
            
             % 1) a random integer to select the object in the list of 
             % coating stl files to be added to the aggregate
             
             if closet.vector_coat_sizes_mu>63
                order_RND_yn = 1;
             else
                order_RND_yn = 2;
             end
                
            
             % 2) scaling of the selected object to the actual size desired by the user 
             singleCoat = realObjectCreator(plusier_STLs.coat(order_RND_yn(1)), ...
                                            closet.vector_coat_sizes_mu(id_particle));
            
             % ----------------- Scaling of the coating: END --------------
            
             parfor id_angle = 1:closet.N_Euler_triplets                   % parallel loop
        
                % Remember the order of the Euler's angles: fi, theta, psi
                vec_rotation = [getRandomNumber(0, 180), ...
                                getRandomNumber(0, 180), ...
                                getRandomNumber(0, 360)];                        
                
               
                % Here we create the ACTUAL PARTICLE object: rotated
                % AND translated
                p_Object = [];
                idxs_blob = [];                                                          

                %%%%%%%%%%%%%%%%%%%%%      LOOP 4     %%%%%%%%%%%%%%%%%%%%%
                
                p_upper_x = start_p.x + closet.vector_coat_sizes_mu(id_particle) * x_dir;
                p_upper_y = start_p.y + closet.vector_coat_sizes_mu(id_particle) * y_dir;
                p_upper_z = start_p.z + closet.vector_coat_sizes_mu(id_particle) * z_dir;
                
                p_Object = createPObject(singleCoat.fv, singleCoat.sphere_struct, ...
                                         closet, vec_rotation, id_particle, ...
                                         [p_upper_x p_upper_y p_upper_z]); 
                
%                 [p_upper_x, p_upper_y, p_upper_z, p_Object] = checkInitialEmptyIntersection(p_Object, s_Blob, closet, singleCoat, vec_rotation, ...
%                                                                            p_upper_x, p_upper_y, p_upper_z, id_particle, ...
%                                                                            x_dir, y_dir, z_dir, start_p);                                                         
                                                                       
                [p_Object, p_upper_x, p_upper_y, p_upper_z] = moveObjectInward(p_Object, s_Blob, closet, ...
                                                                               singleCoat, vec_rotation, ...
                                                                               p_upper_x, p_upper_y, p_upper_z, id_particle, ...
                                                                                x_dir, y_dir, z_dir); 
                                                                                                                                  
                [p_Object] = moveObjectOutward(p_Object, s_Blob, closet, ...
                                               singleCoat, vec_rotation, ...
                                               p_upper_x, p_upper_y, p_upper_z, id_particle, ...
                                               x_dir, y_dir, z_dir);
                                     
                
                % Here we find the center of mass of the new generated
                % object (N.B. for each new generated object)
                [x_CM, y_CM, z_CM] = findCenterOfMass2(p_Object.sphere_struct.x_c, ...
                                                       p_Object.sphere_struct.y_c, ...
                                                       p_Object.sphere_struct.z_c, ... 
                                                       p_Object.sphere_struct.r);
        
                % This vector has all the distances between the center of
                % mass of the central object and the single new object
                dist_Euler(id_angle) = distCalculator(c_Object.CM(1), c_Object.CM(2), c_Object.CM(3), ...
                                                      x_CM, y_CM, z_CM);
                
                % Here we create the prototype object "single"
                % DEBUG: questa funzione puo' essere tolta e attacchi
                % direttamente l'oggetto p_Object alla struttura dati
                % "object_All_Eulers"
                object_Single = createSingleObjectXShot(p_Object, id_angle);
                object_All_Eulers(id_angle).s_Object = object_Single;
                p_Object = [];
        
            end % End of the parallel loop over each rotation of the object
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            %%%%%         END LOOP OVER ALL THE EULER ANGLES          %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                                        % *** %
                                                                
            % Here we find the object closest to the center of the figure
            % within all the Euler angles
            [~, id_dist_min_Euler] = min(dist_Euler);
            
            % Each ray in the cone is related to the minimum Euler angle 
            object_Within_Cone(id_single_ray).s_Object = object_All_Eulers(id_dist_min_Euler).s_Object;
            
            % This is the distance from the center of mass applied ONLY to the closest object along one single ray 
            [x_CM, y_CM, z_CM] = findCenterOfMass2(object_Within_Cone(id_single_ray).s_Object.sphere_struct.x_c, ...
                                                   object_Within_Cone(id_single_ray).s_Object.sphere_struct.y_c, ...
                                                   object_Within_Cone(id_single_ray).s_Object.sphere_struct.z_c, ...
                                                   object_Within_Cone(id_single_ray).s_Object.sphere_struct.r);
                       
            % This is the vector where for each ray you have the distance of the object relative to the center of mass 
            dist_rays(id_single_ray) = distCalculator(c_Object.CM(1), c_Object.CM(2), c_Object.CM(3), x_CM, y_CM, z_CM);                                  
            
         end                                                               % End of the for loop over each ray-tracing
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
         %%%%%    END LOOP OVER ALL THE RAYS FOR A SINGLE PARTICLE    %%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                                % *** %
        
         % In the structure "object_Within_Cone.s_Object" you store all the
         % closest objects for each ray emitted
         [~, id_dist_ray] = min(dist_rays);
        
         g_Object = updateGlobalObject(g_Object, object_Within_Cone(id_dist_ray), id_particle+1);
         s_Blob = mergeObjects(s_Blob, g_Object(id_particle+1).s_Object);
         
         
         
         [output, ~, ~, ~] = calculatePorosity(closet, s_Blob, g_Object, id_particle);
         closet.vettore_porosita(id_particle) = output.porosity;
%         closet.vettore_porosita(id_particle) = 1;
        
        
    end % End of the for loop over each single particle
    toc                 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%          POROSITY CALCULATION USING THE CONVEX HULL         %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [output, C, DT_external, front_points] = calculatePorosity_old(closet, s_Blob, g_Object);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                       PLOTTING PART                         %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     figure(1)
%     
%     for i=1:closet.N_particles+1
%         renderSTL(g_Object(i).s_Object.fv, closet)
%         hold on
%     end
%     axis equal
%     
%     figure(2)
%     
%     for i=1:closet.N_particles+1
%         renderSTL(g_Object(i).s_Object.fv, closet)
%         hold on
%     end
%     axis equal
%     
% %     figure(3)
%     trisurf(C,DT_external.Points(:,1),DT_external.Points(:,2),DT_external.Points(:,3), 'FaceColor','cyan')
%     hold on
%     plot3(front_points(:,1), front_points(:,2), front_points(:,3),'.r', 'MarkerSize', 20)
%     axis equal
%     
%     figure(3)
%     plot3(front_points(:,1), front_points(:,2), front_points(:,3), '.r', 'MarkerSize', 20)
%     hold on
%     plot3(x_CM, y_CM, z_CM, '+r', 'MarkerSize', 30)
%     axis equal
%     
%     figure(4)
%     plotBlob(s_Blob)

%     profile viewer


end

% ----------------------------------------------------------------------- %
% Initialize the coating particles

function [plusier_STLs, closet] = elaborateInitialStructure(cored_struct_original, closet)
   
    % This function is in charge of translating the core to the center of
    % coordinates (x=0, y=0, z=0) and to create the new structure
    % "plusier_STLs", which will be used in the code to not overwrite the
    % input structure. 
    % The structure of plusier_STLs is as follows:
    % - plusier_STLs.core.fv
    % - plusier_STLs.core.sphere_struct
    % - plusier_STLs.core.Feret
    % - plusier_STLs.coat(1).fv
    % - plusier_STLs.coat(1).sphere_struct
    % - plusier_STLs.coat(1).Feret
    % - plusier_STLs.coat(closet.n_obj_coat).fv
    % - plusier_STLs.coat(closet.n_obj_coat).sphere_struct
    % - plusier_STLs.coat(closet.n_obj_coat).Feret
    
                                % *** %
    
    n_obj_tot = length(cored_struct_original);                             % Number of the total objects (stl files) set up by the user to create the aggregate
    closet.n_obj_coat = n_obj_tot-1;                                       % Number of objects dedicated to the coating part  
    
    % Translation of the origin of the coordinates: it is a
    % double step process, since we need to translate both the spherical
    % representation of the object and the stl files associated with it.
    [plusier_STLs.core.fv, plusier_STLs.core.sphere_struct] = center2Zero(cored_struct_original(1).fv, cored_struct_original(1).sphere_struct);
    
    plusier_STLs.core.Feret = getFeret(plusier_STLs.core.fv.vertices, closet);
  
    for i=1:closet.n_obj_coat
        
        [plusier_STLs.coat(i).fv, plusier_STLs.coat(i).sphere_struct] = center2Zero(cored_struct_original(i+1).fv, cored_struct_original(i+1).sphere_struct);
        plusier_STLs.coat(i).Feret = getFeret(plusier_STLs.coat(i).fv.vertices, closet);
        
    end
    
end

function [closet] = initializeEmptyMatrices(closet)

    closet.dist_Euler = zeros(closet.N_Euler_triplets,1);
    closet.object_All_Eulers(closet.N_Euler_triplets).s_Object = [];
    closet.object_Within_Cone(closet.N_raysXSA).s_Object = [];
    closet.dist_rays = zeros(closet.N_raysXSA,1);
    closet.g_Object(closet.N_particles+1).s_Object = [];
    closet.memory_start_p = zeros(3,1);
    closet.CM_prov = zeros(1,3);
    translation_trigger = zeros(closet.N_Euler_triplets,1);
    translation_counter = ones(closet.N_Euler_triplets,1);
    n_of_intersections = zeros(closet.N_Euler_triplets,1);

end

function welcomeMessage
    
    clc
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('                           WELCOME to SCARLET v1.0! \n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('SCARLET is a platform to create virtual aggregates of arbitrary shapes.\n')
    fprintf('Enjoy SCARLET for science or ... fun! ;) ')
    fprintf('\n')
    fprintf('\n')
    fprintf('Authors: Eduardo Rossi, Costanza Bonadonna');
    fprintf('\n')
    fprintf('2020 - All rights are reserved')
    pause(3)
    clc

end

% ----------------------------------------------------------------------- %
% Create virtual objects

function central_Object = createCObject(fv, sphere_struct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% This function creates the central object %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    central_Object.fv = fv;
    central_Object.sphere_struct = sphere_struct;
    [x_CM, y_CM, z_CM] = findCenterOfMass2(central_Object.sphere_struct.x_c, central_Object.sphere_struct.y_c, central_Object.sphere_struct.z_c, central_Object.sphere_struct.r);
    central_Object.CM(1) = x_CM;
    central_Object.CM(2) = y_CM;
    central_Object.CM(3) = z_CM;
    central_Object.N_spheres = length(sphere_struct.id);
    
end

function particle_Object = createPObject(fv, sphere_struct, closet, vec_rotation, id_object, start_p)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Create the single object relative to a particle. All the fields
    %%%%% will be updated during the evolution of the object, but here you
    %%%%% simply create the object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Ps. vec_rotation --> [fi, theta, psi]
    
    particle_Object.id_p = id_object;
    particle_Object.sphere_struct = sphere_struct;
    
    particle_Object.fv.faces = fv.faces;
    particle_Object.fv.vertices(:,1) = fv.vertices(:,1);
    particle_Object.fv.vertices(:,2) = fv.vertices(:,2);
    particle_Object.fv.vertices(:,3) = fv.vertices(:,3);
    
%     % Scaling of the object
%     particle_Object = scaleObjectFunction(particle_Object, closet.vector_scaling(id_object));
    
    % Rotation of the object
    particle_Object = objectRotated(particle_Object, vec_rotation(1), vec_rotation(2), vec_rotation(3));
    
    % Translation of the object
    particle_Object.sphere_struct.x_c = particle_Object.sphere_struct.x_c' + start_p(1);
    particle_Object.sphere_struct.y_c = particle_Object.sphere_struct.y_c' + start_p(2);
    particle_Object.sphere_struct.z_c = particle_Object.sphere_struct.z_c' + start_p(3);
    
    particle_Object.fv.vertices(:,1) = particle_Object.fv.vertices(:,1) + start_p(1);
    particle_Object.fv.vertices(:,2) = particle_Object.fv.vertices(:,2) + start_p(2);
    particle_Object.fv.vertices(:,3) = particle_Object.fv.vertices(:,3) + start_p(3);
    
    % Center of mass of the object
    
    [x_CM, y_CM, z_CM] = findCenterOfMass2(sphere_struct.x_c, sphere_struct.y_c, sphere_struct.z_c, sphere_struct.r);
    
    particle_Object.CM(1) = x_CM;
    particle_Object.CM(2) = y_CM;
    particle_Object.CM(3) = z_CM;
    
   
end

function object_Single = createSingleObjectXShot(particle_Object, id_shot)

    object_Single.sphere_struct = particle_Object.sphere_struct;
    object_Single.fv.vertices = particle_Object.fv.vertices;
    object_Single.fv.faces = particle_Object.fv.faces;
    object_Single.id = id_shot;
    object_Single.CM = particle_Object.CM;
    
end

function s_Blob = initializeBlob(central_Object)

    % Initialization of the structure that contains all the spheres of all
    % the aggregates objects

    s_Blob.sphere_struct.r = central_Object.sphere_struct.r;
    s_Blob.sphere_struct.x_c = central_Object.sphere_struct.x_c;
    s_Blob.sphere_struct.y_c = central_Object.sphere_struct.y_c;
    s_Blob.sphere_struct.z_c = central_Object.sphere_struct.z_c;
    
    [s_Blob.x_CM, s_Blob.y_CM, s_Blob.z_CM] = findCenterOfMass2(central_Object.sphere_struct.x_c, ...
                                                                central_Object.sphere_struct.y_c, ...
                                                                central_Object.sphere_struct.z_c, ...
                                                                central_Object.sphere_struct.r);
                                                            
   s_Blob.n_tot_spheres = length(central_Object.sphere_struct.r);
    
end

function s_Blob = mergeObjects(global_Object, single_Object)

    n_old_spheres = length(global_Object.sphere_struct.r);
    n_new_spheres = length(single_Object.sphere_struct.r);
    n_tot_spheres = n_old_spheres+n_new_spheres;
    
    s_Blob.sphere_struct.r(n_tot_spheres) = -1;
    s_Blob.sphere_struct.x_c(n_tot_spheres) = -1;
    s_Blob.sphere_struct.y_c(n_tot_spheres) = -1;
    s_Blob.sphere_struct.z_c(n_tot_spheres) = -1;
    
    % Placing the old spheres
    s_Blob.sphere_struct.r(1:n_old_spheres) = global_Object.sphere_struct.r;
    s_Blob.sphere_struct.x_c(1:n_old_spheres) = global_Object.sphere_struct.x_c;
    s_Blob.sphere_struct.y_c(1:n_old_spheres) = global_Object.sphere_struct.y_c;
    s_Blob.sphere_struct.z_c(1:n_old_spheres) = global_Object.sphere_struct.z_c;
    
    % Placing the new spheres
    s_Blob.sphere_struct.r(n_old_spheres+1:end) = single_Object.sphere_struct.r;
    s_Blob.sphere_struct.x_c(n_old_spheres+1:end) = single_Object.sphere_struct.x_c;
    s_Blob.sphere_struct.y_c(n_old_spheres+1:end) = single_Object.sphere_struct.y_c;
    s_Blob.sphere_struct.z_c(n_old_spheres+1:end) = single_Object.sphere_struct.z_c;
    
    % Center of mass of the new global object
    [s_Blob.x_CM, s_Blob.y_CM, s_Blob.z_CM] = findCenterOfMass2(s_Blob.sphere_struct.x_c, ...
                                                                s_Blob.sphere_struct.y_c, ...
                                                                s_Blob.sphere_struct.z_c, ...
                                                                s_Blob.sphere_struct.r);
    s_Blob.n_tot_spheres = n_tot_spheres;

end

function global_Object = updateGlobalObject(global_Object, particle_Object_single, idx)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%         This function updates the global object             %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In practice: it merges the central object with all the final object %
    % selected within each solid angle.                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global_Object(idx).s_Object = particle_Object_single.s_Object;

end


% ----------------------------------------------------------------------- %
% Rotation, translation, scaling functions

function [fv, sphere_struct] = center2Zero(fv_original, sphere_struct_original)

    % Function to translate the coordinates of the object to (0,0,0)

    % Center of mass of the object
    [x_CM, y_CM, z_CM] = findCenterOfMass2(sphere_struct_original.x_c, ...
                                           sphere_struct_original.y_c, ...
                                           sphere_struct_original.z_c, ...
                                           sphere_struct_original.r);
       
    % Translation of the points contained in the spherical representation of the object 
    [sphere_struct.x_c, sphere_struct.y_c, sphere_struct.z_c] = initialTranslation(x_CM, y_CM, z_CM, ...
                                                                                   sphere_struct_original.x_c, ...
                                                                                   sphere_struct_original.y_c, ...
                                                                                   sphere_struct_original.z_c);
    % Creation of the new "sphere_struct" translated 
    sphere_struct.r = sphere_struct_original.r;
    sphere_struct.vol = sphere_struct_original.vol;
    sphere_struct.id = sphere_struct_original.id;
    
    % Translation of the points contained in the vertices of the stl file
    [new_vertices_x, new_vertices_y, new_vertices_z] = initialTranslation(x_CM, y_CM, z_CM, ...
                                                                          fv_original.vertices(:,1), ...
                                                                          fv_original.vertices(:,2), ...
                                                                          fv_original.vertices(:,3));
    tab(:,1) = new_vertices_x;
    tab(:,2) = new_vertices_y;
    tab(:,3) = new_vertices_z;
    
    % creation of the translated version of the stl file
    fv.vertices = tab;
    fv.faces = fv_original.faces;
    

end

function [particle_Object] = scaleObjectFunction(particle_Object_not_scaled, c_factor)

    particle_Object.fv.faces = particle_Object_not_scaled.fv.faces;
    particle_Object.fv.vertices = c_factor * particle_Object_not_scaled.fv.vertices;
    
    particle_Object.sphere_struct.x_c = c_factor * particle_Object_not_scaled.sphere_struct.x_c;
    particle_Object.sphere_struct.y_c = c_factor * particle_Object_not_scaled.sphere_struct.y_c;
    particle_Object.sphere_struct.z_c = c_factor * particle_Object_not_scaled.sphere_struct.z_c;
    particle_Object.sphere_struct.r = c_factor * particle_Object_not_scaled.sphere_struct.r;
    particle_Object.sphere_struct.id = particle_Object_not_scaled.sphere_struct.id;
    
    particle_Object.Feret = particle_Object_not_scaled.Feret;

end

function object = realObjectCreator(object, d_real_micron)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% This function creates a real scaled object according to the Feret
    %%%%% diameter chosen in the input file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c_factor = d_real_micron/object.Feret;
    [object] = scaleObjectFunction(object, c_factor);
    
end

function [particle_Object] = objectRotated(particle_Object, fi, theta, psi)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function rotates an object "fv" in the space according to the
    % three Euler angles fi, theta and psi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculation of the Euler matrix (common to the FV object and the
    % spheres
    Euler_matrix = EulerRotationMatrix(deg2rad(fi), deg2rad(theta), deg2rad(psi));
    
    %%%%%%%%%%%%%%%%%%%%%% Rotation of the FV object %%%%%%%%%%%%%%%%%%%%%%
    
    % NB: points_matrix must have a size 3 X number of points
    points_matrix = particle_Object.fv.vertices';
    
    % New coordinates, rotated
    particle_Object.fv.vertices = (Euler_matrix*points_matrix)';
    
    
    %%%%%%%%%%%%%%%%%%%%%%% Rotation of the spheres %%%%%%%%%%%%%%%%%%%%%%%
    
    points_matrix = [];
    
    % NB: points_matrix must have a size 3 X number of points
    points_matrix(1,:) = particle_Object.sphere_struct.x_c;
    points_matrix(2,:) = particle_Object.sphere_struct.y_c;
    points_matrix(3,:) = particle_Object.sphere_struct.z_c;
    
    % New coordinates, rotated
    points_rotated = (Euler_matrix*points_matrix)';
    
    % New structure
    particle_Object.sphere_struct.x_c = points_rotated(:,1);
    particle_Object.sphere_struct.y_c = points_rotated(:,2);
    particle_Object.sphere_struct.z_c = points_rotated(:,3);

end

function [sphere_struct_out] = spheresRotated(sphere_struct, fi, theta, psi)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function rotates all the spheres that contribute to create the
    % skeleton of the fv file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % NB: points_matrix must have a size 3 X number of points
    points_matrix(1,:) = sphere_struct.x_c;
    points_matrix(2,:) = sphere_struct.y_c;
    points_matrix(3,:) = sphere_struct.z_c;

    % Calculation of the Euler matrix
    Euler_matrix = EulerRotationMatrix(deg2rad(fi), deg2rad(theta), deg2rad(psi));
    
    % New coordinates, rotated
    points_rotated = (Euler_matrix*points_matrix)';
    
    % New structure
    sphere_struct_out.x_c = points_rotated(:,1);
    sphere_struct_out.y_c = points_rotated(:,2);
    sphere_struct_out.z_c = points_rotated(:,3);
    sphere_struct_out.r = sphere_struct.r;
    sphere_struct_out.id = sphere_struct.id;

end

function Euler_matrix = EulerRotationMatrix(fi, theta, psi)

   A_11 = cos(psi).*cos(fi) - cos(theta) .* sin(fi) .* sin(psi);
   A_12 = cos(psi).*sin(fi) + cos(theta) .* cos(fi) .* sin(psi);
   A_13 = sin(psi).*sin(theta);
   
   A_21 = -sin(psi).*cos(fi) - cos(theta) .* sin(fi) .* cos(psi);
   A_22 = -sin(psi).*sin(fi) + cos(theta) .* cos(fi) .* cos(psi);
   A_23 = cos(psi).*sin(theta);
   
   A_31 = sin(theta).*sin(fi);
   A_32 = -sin(theta).*cos(fi);
   A_33 = cos(theta);
   
   Euler_matrix = [A_11 A_12 A_13; A_21 A_22 A_23; A_31 A_32 A_33];
   
end

function [new_vec_x, new_vec_y, new_vec_z] = initialTranslation(x_cm, y_cm, z_cm, vec_x, vec_y, vec_z)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Translation of the object in order to overlap the CoM with (0,0,0) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end_point_x = 0;
    end_point_y = 0;
    end_point_z = 0;
    
    % Vector translation
    dx = end_point_x-x_cm;
    dy = end_point_y-y_cm;
    dz = end_point_z-z_cm;
    
    % Translated points
    new_vec_x = vec_x+dx;
    new_vec_y = vec_y+dy;
    new_vec_z = vec_z+dz;
    
%     % NOT TRANLSATED OPTION
%     new_vec_x = vec_x;
%     new_vec_y = vec_y;
%     new_vec_z = vec_z;
    
end

function Feret = getFeret(vertices_table, closet)

    % Goal of the function: to measure the maximum distance between the
    % points of the object. For "object" we mean one single stl file or
    % more.
    
    % Input of the function: "vertices_table", which is a table N_rowsX3
    % containing all the vertices of the objects, i.e. the points.
    
    % Additional description: in order to speed-up the "pdist" function we
    % take only a random subset of the total points present in the object.
    % The number of the points is passed from outside by means of "closet"
    
    dim = size(vertices_table);
    idxs_rnd = randperm(dim(1));
    idxs_rnd_net = idxs_rnd(1:closet.n_points4Feret);
    tab_4_Feret = vertices_table(idxs_rnd_net,:);
    Feret = max(pdist(tab_4_Feret));
    
end

function [p_Object, p_upper_x, p_upper_y, p_upper_z] = moveObjectInward(p_Object, s_Blob, closet, singleCoat, vec_rotation, ...
                                                                        p_upper_x, p_upper_y, p_upper_z, id_particle, ...
                                                                        x_dir, y_dir, z_dir)
   
        translation_trigger = 0;
        translation_counter = 1;

        while (translation_trigger < 1)
                                                  
        n_of_intersections = countIntersections2020(p_Object, s_Blob);
                        
        if n_of_intersections==0 
                            
            p_upper_x = p_upper_x - closet.vector_coat_sizes_mu(id_particle) * x_dir*translation_counter * closet.delta;
            p_upper_y = p_upper_y - closet.vector_coat_sizes_mu(id_particle) * y_dir*translation_counter * closet.delta;
            p_upper_z = p_upper_z - closet.vector_coat_sizes_mu(id_particle) * z_dir*translation_counter * closet.delta;
                            
            p_Object = createPObject(singleCoat.fv, singleCoat.sphere_struct, closet, vec_rotation, id_particle, ...
                                     [p_upper_x p_upper_y p_upper_z]); 
                            
        else
            translation_trigger = 1;
%             translation_trigger(id_angle) = 1;
        end
                                                                     
        if (translation_counter > closet.translation_iter_max)
            disp('WARNING: reached MAX multiple iterations in the while loop within SCARLETAggrAggr4General')
            translation_trigger = 1;
%             translation_trigger(id_angle) = 1;
        end
        
        % This integer is increased each while iteration
        translation_counter = translation_counter + 1;
                    
        end

end

function [p_Object] = moveObjectOutward(p_Object, s_Blob, closet, singleCoat, vec_rotation, ...
                                       p_upper_x, p_upper_y, p_upper_z, id_particle, ...
                                       x_dir, y_dir, z_dir)
   
        translation_trigger = 0;
        translation_counter = 1;

        while (translation_trigger < 1)
                                                  
        n_of_intersections = countIntersections2020(p_Object, s_Blob);
                        
        if n_of_intersections>0 
                            
            p_upper_x = p_upper_x + closet.vector_coat_sizes_mu(id_particle) * x_dir*translation_counter * closet.delta_2;
            p_upper_y = p_upper_y + closet.vector_coat_sizes_mu(id_particle) * y_dir*translation_counter * closet.delta_2;
            p_upper_z = p_upper_z + closet.vector_coat_sizes_mu(id_particle) * z_dir*translation_counter * closet.delta_2;
                            
            
            p_Object = createPObject(singleCoat.fv, singleCoat.sphere_struct, closet, vec_rotation, id_particle, ...
                                     [p_upper_x p_upper_y p_upper_z]); 
                            
        else
            translation_trigger = 1;
%             translation_trigger(id_angle) = 1;
        end
                                                                     
        if (translation_counter > closet.translation_iter_max)
            disp('WARNING: reached MAX multiple iterations in the while loop within SCARLETAggrAggr4General')
            translation_trigger = 1;
%             translation_trigger(id_angle) = 1;
        end
        
        % This integer is increased each while iteration
        translation_counter = translation_counter + 1;
                    
        end

end

function [p_upper_x, p_upper_y, p_upper_z, p_Object] = checkInitialEmptyIntersection(p_Object, s_Blob, closet, singleCoat, vec_rotation, ...
                                                                                     p_upper_x, p_upper_y, p_upper_z, id_particle, ...
                                                                                     x_dir, y_dir, z_dir, start_p)
                                                                       
    % This function checks if the initial intersection between the particle
    % and the existing structure is zero. If not, the code moves outward
    % the coating particle 
                                   
    n_of_intersections_check = 1;
    inner_counter = 1;
    
    while (n_of_intersections_check>0)
                                        
        n_of_intersections_check = countIntersections2020(p_Object, s_Blob);
                    
        if n_of_intersections_check>0
            p_upper_x = start_p.x + inner_counter*closet.vector_coat_sizes_mu(id_particle) * x_dir;
            p_upper_y = start_p.y + inner_counter*closet.vector_coat_sizes_mu(id_particle) * y_dir;
            p_upper_z = start_p.z + inner_counter*closet.vector_coat_sizes_mu(id_particle) * z_dir;
            
            p_Object = createPObject(singleCoat.fv, singleCoat.sphere_struct, ...
                                 closet, vec_rotation, id_particle, ...
                                 [p_upper_x p_upper_y p_upper_z]);  
            
            
            inner_counter = inner_counter+1;
        end
                    
     end

end

% ----------------------------------------------------------------------- %
% Interesection routines

function dist = distCalculator(x1, y1, z1, x2, y2, z2)

    dist = sqrt( (x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2 );

end

function a = intersection2SpheresAnalytic(R, r, d)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Analytical intersection of two spheres of equations: 
    %%% - R is the radius of the old sphere (the one belonging to the one
    %%%  already placed in the location)
    %%% - r is the radis of the new sphere (the one belonging to the object
    %%% to be placed)
    %%% - d distance between the centers of the two spheres
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Radius of the intersecting circonference
    a = 1./(2.*d) .* sqrt( 4 .* d.^2 .* R.^2 - (d.^2 - r.^2 + R.^2).^2 ); 
    
end

function n_of_intersections = countIntersections2020(particle_Object_new, particle_Object_old)

    % We compare the "new" object (its internal spheres) with those already
    % existing
    
    N_new_spheres = length(particle_Object_new.sphere_struct.r);
    N_old_spheres = length(particle_Object_old.sphere_struct.r);
    
    vector_intersections = zeros(N_new_spheres, 1);
    r = zeros(1,1);
    R = r;
    d_x = zeros(1,N_old_spheres);
    d_y = d_x;
    d_z = d_x;
    d = d_x;
    vector = d_x; 
    
    % Here you compare one single sphere of the new object with all the
    % other spheres already placed
    %tic
    for i=1:N_new_spheres
        
        r = particle_Object_new.sphere_struct.r(i);
        R = particle_Object_old.sphere_struct.r;
        
        d_x = particle_Object_old.sphere_struct.x_c - particle_Object_new.sphere_struct.x_c(i);
        d_y = particle_Object_old.sphere_struct.y_c - particle_Object_new.sphere_struct.y_c(i);
        d_z = particle_Object_old.sphere_struct.z_c - particle_Object_new.sphere_struct.z_c(i);
    
        d = sqrt(d_x.^2 + d_y.^2 + d_z.^2);
        
        vector = intersection2SpheresAnalytic(R, r, d);
        idxs = sum(imag(vector)==0);
        
        vector_intersections(i) = sum(imag(vector) == 0);
        
    end
    %toc
   
    n_of_intersections = sum(vector_intersections);
    
end

function [p_inter_max, p_inter_min] = intersectionLineSphere(p_Object, L_px, L_vx, L_py, L_vy, L_pz, L_vz)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Intersection of one line with multiple spheres (ray-tracing intersection)
    %%%%% p_inter_max = it is the most "positive" intersecting point,
    %%%%% therefore the most "exterior one"
    %%%%%
    %%%%% p_inter_min = it is the most "negative" intersecting point,
    %%%%% therefore the most "opposite one"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S_cx = p_Object.sphere_struct.x_c;
    S_cy = p_Object.sphere_struct.y_c;
    S_cz = p_Object.sphere_struct.z_c;
    S_r = p_Object.sphere_struct.r;
    
    A = (L_vx .* L_vx + L_vy .* L_vy + L_vz .* L_vz);
    
    B = 2 .* (L_px .* L_vx + L_py .* L_vy + L_pz .* L_vz - L_vx .* S_cx - ...
              L_vy .* S_cy - L_vz .* S_cz);
    
    C = L_px .* L_px - 2 .* L_px .* S_cx + S_cx .* S_cx + L_py .* L_py - ...
        2 .* L_py .* S_cy + S_cy .* S_cy + L_pz .* L_pz - 2 .* L_pz .* S_cz  + ...
        S_cz .* S_cz - S_r .* S_r;
    
    D = B .* B - 4 .* A .* C;
    
    % First set of solutions
    t_1 = (-B - sqrt(D)) ./ (2 .* A);
    
    % Second set of solutions
    t_2 = (-B + sqrt(D)) ./ (2 .* A);
    
    % First set of solutions
    sol_x_1 = L_px + t_1.*L_vx;
    sol_y_1 = L_py + t_1.*L_vy;
    sol_z_1 = L_pz + t_1.*L_vz;
    
    p_1_x_imag = imag(sol_x_1);
    p_1_y_imag = imag(sol_y_1);
    p_1_z_imag = imag(sol_z_1);
    
    p_1_sum = p_1_x_imag + p_1_y_imag + p_1_z_imag;
    
    % "idxs_1" represents the vector with indexes of real points
    idxs_1 = find(p_1_sum==0);
    
    p_inter.x = sol_x_1(idxs_1);
    p_inter.y = sol_y_1(idxs_1);
    p_inter.z = sol_z_1(idxs_1);
    
    % Second set of solutions
    sol_x_2 = L_px + t_2.*L_vx;
    sol_y_2 = L_py + t_2.*L_vy;
    sol_z_2 = L_pz + t_2.*L_vz;
    
    p_2_x_imag = imag(sol_x_2);
    p_2_y_imag = imag(sol_y_2);
    p_2_z_imag = imag(sol_z_2);
    
    p_2_sum = p_2_x_imag+p_2_y_imag+p_2_z_imag;
    
    % "idxs_1" represents the vector with indexes of real points
    idxs_2 = find(p_2_sum==0);
    
    p_inter.x = [p_inter.x sol_x_2(idxs_2)];
    p_inter.y = [p_inter.y sol_y_2(idxs_2)];
    p_inter.z = [p_inter.z sol_z_2(idxs_2)];
    dot_products = zeros(length(p_inter.x),1);
    
    for i=1:length(p_inter.x)
        
        dot_products(i) = dot([L_vx L_vy L_vz], [p_inter.x(i); p_inter.y(i); p_inter.z(i)]);
        
    end
    
    [~, id_max] = max(dot_products);
    [~, id_min] = min(dot_products);
    
    p_inter_max.x = p_inter.x(id_max);
    p_inter_max.y = p_inter.y(id_max);
    p_inter_max.z = p_inter.z(id_max);
    
    p_inter_min.x = p_inter.x(id_min);
    p_inter_min.y = p_inter.y(id_min);
    p_inter_min.z = p_inter.z(id_min);

end

function id = IDdistCalculatorVector(x_p, y_p, z_p, x_vec, y_vec, z_vec)

    dist = sqrt( (x_p-x_vec).^2 + (y_p-y_vec).^2 + (z_p-z_vec).^2 );
    [min_val, id] = min(dist); 

end

function [front_points, x_CM, y_CM, z_CM]  = findPoints4ExternalSurface(s_Blob)
    
    n_points = 10000;
    versor_table = generateUnitVector(n_points);
    cone_central_dir = zeros(3,1);
    front_points = zeros(size(versor_table));
    points.x = 0;
    points.y = 0;
    points.z = 0;
    
    [x_CM, y_CM, z_CM] = findCenterOfMass2(s_Blob.sphere_struct.x_c, s_Blob.sphere_struct.y_c, s_Blob.sphere_struct.z_c, s_Blob.sphere_struct.r);
          
    for i=1:n_points
        
        cone_central_dir = [versor_table(i,1); versor_table(i,2); versor_table(i,3)];
        
        [points, ~] = intersectionLineSphere(s_Blob, ...
                                             x_CM, cone_central_dir(1), ...
                                             y_CM, cone_central_dir(2), ...
                                             z_CM, cone_central_dir(3));
        
        if isempty(points.x)>0
            
            id = IDdistCalculatorVector(x_CM, y_CM, z_CM, ...
                                        s_Blob.sphere_struct.x_c, s_Blob.sphere_struct.y_c, s_Blob.sphere_struct.z_c);
            
            [points, ~] = intersectionLineSphere(s_Blob, ...
                                                 s_Blob.sphere_struct.x_c(id), cone_central_dir(1), ...
                                                 s_Blob.sphere_struct.y_c(id), cone_central_dir(2), ...
                                                 s_Blob.sphere_struct.z_c(id), cone_central_dir(3));
            
        end
        
               
        front_points(i,1) = points.x;
        front_points(i,2) = points.y;
        front_points(i,3) = points.z;
   
    end

    
end
% ----------------------------------------------------------------------- %
% Plotting routines

function plotSpheres(particleObject, closet)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% This function plots the spheres %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     
     for i = 1:length(particleObject.sphere_struct.r)
        
        [points_x_c, points_y_c, points_z_c] = singleSphere4Plot(particleObject.sphere_struct.r(i), ...
                                                                 particleObject.sphere_struct.x_c(i), ...
                                                                 particleObject.sphere_struct.y_c(i), ...
                                                                 particleObject.sphere_struct.z_c(i));
        
        hl = surf(points_x_c, points_y_c, points_z_c);
        colormap winter
        shading interp
        alpha(hl, closet.alpha_transparency)
        axis equal
        view([closet.view_angle_1 closet.view_angle_2]);
        xlabel('x axis')
        ylabel('y axis')
        zlabel('z axis')
        drawnow
        hold on
        
     end

%      hold off
     
end

function renderSTL(fv, closet)

%     patch(fv,'FaceColor', [0.8 0.8 1.0], ...
%              'EdgeColor', 'none', ...
%              'FaceLighting', 'gouraud', ...
%              'AmbientStrength', 0.85);

    patch(fv,'FaceColor', [1 0 0], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud', ...
          'AmbientStrength', 0.85);  
          
    % Add a camera light, and tone down the specular highlighting
    
%     camlight('headlight');
    lightangle(180,60)
    material metal 
    alpha(closet.alpha_transparency)                                                             % Flag: ho aggiunto questo solo per vedere le sfere dentro
    
    view([closet.view_angle_1 closet.view_angle_2]);
%     view([-141 36]);
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    axis equal
    
end

function renderSTL2(fv)

    patch(fv,'FaceColor', [0.8 0.8 1.0], ...
             'EdgeColor', 'none', ...
             'FaceLighting', 'gouraud', ...
             'AmbientStrength', 0.85);

    
    % Add a camera light, and tone down the specular highlighting
    
%     camlight('headlight');
    lightangle(180,60)
    material metal 
    alpha(0.1)                                                             % Flag: ho aggiunto questo solo per vedere le sfere dentro
    
    view([35 42]);
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    axis equal
    
end

function [points_x_c, points_y_c, points_z_c] = singleSphere4Plot(r, x_c, y_c, z_c)

    % General sphere generator for coating
    [x, y, z] = sphere(25);

    points_x_c = x * r + x_c;
    points_y_c = y * r + y_c;
    points_z_c = z * r + z_c;
    
end

% ----------------------------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mathematical tools subroutines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = deg2rad(x)
    y = x.*pi./180;
end

function y = rad2deg(x)
    y = x.*180./pi;
end

function [x_CM, y_CM, z_CM] = findBarycenter(x_vector, y_vector, z_vector)

    x_CM = mean(x_vector);
    y_CM = mean(y_vector);
    z_CM = mean(z_vector);
    
end

function [x_CM, y_CM, z_CM] = findCenterOfMass2(x_vector, y_vector, z_vector, r)

    m = 4/3*pi*r.^3;
    M_tot = sum(m);
    x_CM = sum(m.*x_vector)/M_tot;
    y_CM = sum(m.*y_vector)/M_tot;
    z_CM = sum(m.*z_vector)/M_tot;
    
end

function num_rnd = getRandomNumber(min_value, max_value)

    num_rnd = min_value + rand(1,1)*(max_value-min_value);

end

function num_rnd = getRandomNumberVec(min_value, max_value, N)

    num_rnd = min_value + rand(1,N)*(max_value-min_value);

end

function output = solidAngleRndPoints(coneAngleDegree, coneDir, N, RNG)

    if ~exist('RNG', 'var') || isempty(RNG)
        RNG = RandStream.getGlobalStream();
    end

    coneAngle = deg2rad(coneAngleDegree);

     % Generate points on the spherical cap around the north pole [1].
    % [1] See https://math.stackexchange.com/a/205589/81266
    z = RNG.rand(1, N) * (1 - cos(coneAngle)) + cos(coneAngle);
    phi = RNG.rand(1, N) * 2 * pi;
    x = sqrt(1-z.^2).*cos(phi);
    y = sqrt(1-z.^2).*sin(phi);

    % If the spherical cap is centered around the north pole, we're done.
    if all(coneDir(:) == [0;0;1])
        output = [x; y; z];
        return;
    end

    % Find the rotation axis `u` and rotation angle `rot` [1]
    u = normc(cross([0;0;1], normc(coneDir)));
    rot = acos(dot(normc(coneDir), [0;0;1]));

    % Convert rotation axis and angle to 3x3 rotation matrix [2]
    % [2] See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    crossMatrix = @(x,y,z) [0 -z y; z 0 -x; -y x 0];
    R = cos(rot) * eye(3) + sin(rot) * crossMatrix(u(1), u(2), u(3)) + (1-cos(rot))*(u * u');

    % Rotate [x; y; z] from north pole to `coneDir`.
    output = R * [x; y; z];

end

function y = normc(x)
    y = bsxfun(@rdivide, x, sqrt(sum(x.^2)));
end

function v = generateUnitVector(n)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%             This function generates n unit vectors          %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    v = randn(n,3);
    v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));

end

function [totalVolume, totalArea] = stlVolumeModified(p,t)

    % Given a surface triangulation, compute the volume enclosed using
    % divergence theorem.
    % Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
    % Input: p: (3xnPoints), t: (3xnTriangles)
    % Output: total volume enclosed, and total area of surface  
    % Author: K. Suresh; suresh@engr.wisc.edu
    % To be used as: [totalVolume, totalArea] = stlVolume(fv.vertices', fv.faces');
    % Compute the vectors d13 and d12
    
    d13 = [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
    d12 = [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
    
    cr = cross(d13, d12, 1);                                               % Cross-product (vectorized)
    
    area = 0.5*sqrt(cr(1,:).^2 + cr(2,:).^2 + cr(3,:).^2);                 % Area of each triangle
    
    totalArea = sum(area);
    
    crNorm = sqrt(cr(1,:).^2 + cr(2,:).^2 + cr(3,:).^2);
    
    zMean = (p(3,t(1,:)) + p(3,t(2,:)) + p(3,t(3,:)))/3;
    
    nz = -cr(3,:)./crNorm;                                                 % z component of normal for each triangle
    
    volume = area.*zMean.*nz;                                              % Contribution of each triangle
    
    totalVolume = sum(volume);                                             % Divergence theorem
    
end

function [output, C, DT_external, front_points] = calculatePorosity_old(closet, s_Blob, g_Object)

    % Porosity calculation using the convex hull volume
    
    [front_points, x_CM, y_CM, z_CM] = findPoints4ExternalSurface(s_Blob);
    
    % Convex hull over the external object
    DT_external = delaunayTriangulation(front_points);
    [C, vol_ext] = convexHull(DT_external);
    
    for i = 1:closet.N_particles+1
        
        [iesim_Volume(i), ~] = stlVolumeModified(g_Object(i).s_Object.fv.vertices', g_Object(i).s_Object.fv.faces');
        
    end
    
    output.vol_int = sum(iesim_Volume);
    output.vol_ext = vol_ext;
    output.porosity = 1-output.vol_int/output.vol_ext;
    output.density = closet.core_density * (1-output.porosity);
    output.size = (6*output.vol_ext/pi)^(1/3);
    output.mass = closet.core_density*output.vol_int;

end

function [output, C, DT_external, front_points] = calculatePorosity(closet, s_Blob, g_Object, n_particles_loop)

    % Porosity calculation using the convex hull volume
    
    [front_points, x_CM, y_CM, z_CM] = findPoints4ExternalSurface(s_Blob);
   
    % Convex hull over the external object
    DT_external = delaunayTriangulation(front_points);
    [C, vol_ext] = convexHull(DT_external);
    
    for i = 1:n_particles_loop+1
%     for i = 1:closet.N_particles+1
        
        [iesim_Volume(i), ~] = stlVolumeModified(g_Object(i).s_Object.fv.vertices', g_Object(i).s_Object.fv.faces');
        
    end
    
    output.vol_int = sum(iesim_Volume);
    output.vol_ext = vol_ext;
    output.porosity = 1-output.vol_int/output.vol_ext;

end
