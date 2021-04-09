function [spheres_struct, fv, datalog] = fromStlToSpheres(stl_file, N_spheres, N_iter_max_per_sphere, barycenterMode, light_mode)

    % Preliminary file to fill the stl file with spheres. This file uses
    % the following subroutines: 
    % - "inpolyhedron.m"
    % - "stlread.m"
   
    tic
    
    % Step 1: importing the STL file for reading the external topography of
    % the object
    
    close all
    
    % Number of total spheres
    data.N_spheres_in = N_spheres;
    
    % Boolean variable to discard inner spheres
    data.light_mode = light_mode;
    
    % Number of iterations to find the i-esime sphere
    data.N_check_2_max = N_iter_max_per_sphere;
    
    fv_not_scaled = stlread(stl_file);                                            
   
    disp('Read the stl file')
    
    dim = size(fv_not_scaled.faces);
    
    figure (1)
    renderSTL(fv_not_scaled);
    axis equal
    
    % Step 2: creation of the data structure. 
    data.N_points_net = dim(1);
    
    % Step3: SCALING
    [data, fv] = STLPoints(data, fv_not_scaled, 1);
    data.fv = fv;
    
    figure (2)
    renderSTL(fv);
    axis equal
    view([0 90])
   
    % Step4: find the barycenter
    [data.x_bary, data.y_bary, data.z_bary] = findCentroid(data.x_net, data.y_net, data.z_net);
        
    % Step5: loop over N (net) random points within the object
    
    N_check_2 = 1;
    condition = 1;
    condition_2 = 1;
    vol_spheres_tot = 0;
    const_1 = 4/3*pi; 
    alpha_transparency = 0.9;
    light_sign = 0;
    
    hold on
    
    if (barycenterMode>0)
      
        disp('Barycenter is used as first sphere')  
        
        % Creation of the first sphere
                
        [radius, pos_min, check] = findClosestPoint(data.x_net, data.y_net, data.z_net, data.x_bary, data.y_bary, data.z_bary, 0);

        spheres_struct.r(condition) = radius; 
        spheres_struct.vol(condition) = 4/3*pi*spheres_struct.r(condition)^3;
        spheres_struct.x_c(condition) = data.x_bary;
        spheres_struct.y_c(condition) = data.y_bary;
        spheres_struct.z_c(condition) = data.z_bary;
        spheres_struct.id(condition) = condition; 
                
        % Total volume of the spheres
        vol_spheres_tot = const_1 * spheres_struct.r(condition)^3 + vol_spheres_tot;
                
        [points_x_c, points_y_c, points_z_c] = createSingleSphere(spheres_struct.r(condition), ...
                                                                  spheres_struct.x_c(condition), ...
                                                                  spheres_struct.y_c(condition), ...
                                                                  spheres_struct.z_c(condition));
                                                                
     
        hl = surf(points_x_c, points_y_c, points_z_c);
                    set(hl,'FaceColor',[1 0 0], ...
                    'FaceAlpha',alpha_transparency,'FaceLighting','gouraud','EdgeColor','none')
                    alpha(hl, alpha_transparency)
                    hold on
                     
        drawnow
                                                                                                       
        condition = condition + 1;
        
    else
        disp('Barycenter is NOT used as first sphere')
        
        while (light_sign<1)
            x_rand = (data.max_x-data.min_x).*rand(1,1) + data.min_x;
            y_rand = (data.max_y-data.min_y).*rand(1,1) + data.min_y;
            z_rand = (data.max_z-data.min_z).*rand(1,1) + data.min_z;
            
            IN = inpolyhedron(data.fv, [x_rand, y_rand, z_rand]);
            
            if (IN == 1)
                
                [radius, pos_min, check] = findClosestPoint(data.x_net, data.y_net, data.z_net, x_rand, y_rand, z_rand, 0);
                
                spheres_struct.r(condition) = radius; 
                spheres_struct.vol(condition) = 4/3*pi*spheres_struct.r(condition)^3;
                spheres_struct.x_c(condition) = x_rand;
                spheres_struct.y_c(condition) = y_rand;
                spheres_struct.z_c(condition) = z_rand;
                spheres_struct.id(condition) = condition; 
                
                vol_spheres_tot = const_1 * spheres_struct.r(condition)^3 + vol_spheres_tot;
                [points_x_c, points_y_c, points_z_c] = createSingleSphere(spheres_struct.r(condition), ...
                                                                  spheres_struct.x_c(condition), ...
                                                                  spheres_struct.y_c(condition), ...
                                                                  spheres_struct.z_c(condition));
              
                hl = surf(points_x_c, points_y_c, points_z_c);
                    set(hl,'FaceColor',[1 0 0], ...
                    'FaceAlpha',alpha_transparency,'FaceLighting','gouraud','EdgeColor','none')
                    alpha(hl, alpha_transparency)
                    hold on
                     
                drawnow
                
                condition = condition + 1;
                
                light_sign = 1;
            end
            
        end
        
    end
    
    
     while condition <= data.N_spheres_in
        
        x_rand = (data.max_x-data.min_x).*rand(1,1) + data.min_x;
        y_rand = (data.max_y-data.min_y).*rand(1,1) + data.min_y;
        z_rand = (data.max_z-data.min_z).*rand(1,1) + data.min_z;
        
        IN = inpolyhedron(data.fv, [x_rand, y_rand, z_rand]);
        
        % If here the random point is inside the polygon
        
        if (IN == 1)
                
            % Centers of all the spheres within the volume (dynamic
            % structure: it increases with time)
                
            x_pos_center_spheres = spheres_struct.x_c;
            y_pos_center_spheres = spheres_struct.y_c;
            z_pos_center_spheres = spheres_struct.z_c;
            r_spheres = spheres_struct.r;
               
            [dists_spheres, pos_min, check] = findClosestPoint(x_pos_center_spheres, y_pos_center_spheres, z_pos_center_spheres, x_rand, y_rand, z_rand, r_spheres);
            
            % At this stage we ignore if the point is internal to an
            % existing sphere or not: check!
                
            if (check==0)
    
                [dists_wall, pos_min, check] = findClosestPoint(data.x_net, data.y_net, data.z_net, x_rand, y_rand, z_rand, 0);
                    
                effective_r = min(dists_spheres, dists_wall);
                    
                temporary.r(N_check_2) = effective_r;
                temporary.x_c(N_check_2) = x_rand;
                temporary.y_c(N_check_2) = y_rand;
                temporary.z_c(N_check_2) = z_rand;
                            
                if (N_check_2 > data.N_check_2_max)
                    
                    [dist_net, pos_max] = max(temporary.r);
                        
                    spheres_struct.r(condition) = temporary.r(pos_max); 
                    spheres_struct.vol(condition) = 4/3*pi*spheres_struct.r(condition)^3;
                    spheres_struct.x_c(condition) = temporary.x_c(pos_max); 
                    spheres_struct.y_c(condition) = temporary.y_c(pos_max); 
                    spheres_struct.z_c(condition) = temporary.z_c(pos_max); 
                    spheres_struct.id(condition) = condition; 
                    
                    % Total volume of the spheres
                    vol_spheres_tot = const_1 * spheres_struct.r(condition)^3 + vol_spheres_tot;
                    
                    [points_x_c, points_y_c, points_z_c] = ...
                                                           createSingleSphere(spheres_struct.r(condition), ...
                                                                              spheres_struct.x_c(condition), ...
                                                                              spheres_struct.y_c(condition), ...
                                                                              spheres_struct.z_c(condition));
                    
                    hl = surf(points_x_c, points_y_c, points_z_c);
                    set(hl,'FaceColor',[1 0 0], ...
                    'FaceAlpha',alpha_transparency,'FaceLighting','gouraud','EdgeColor','none')
                    alpha(hl, alpha_transparency)
                    hold on
                     
                    drawnow
                    condition = condition + 1
                    dists = [];
                    check = [];
                    pos_min = [];
                    temporary = [];
                    N_check_2 = 1;
                        
               end
                    
               N_check_2 = N_check_2+1;
                    
           end
        end
        
    end % End While
    
    axis equal
    renderSTL(fv);
    view([0 90])
    
    if data.light_mode>0
        % Avoid the internal spheres
        spheres_struct = filterFinalObject(spheres_struct, fv, N_spheres, alpha_transparency);
    else
        % Not avoud the internal spheres
        spheres_struct = spheres_struct;
    end

   
     % STL volume calculation
    [totalVolume, totalArea] = stlVolume(fv.vertices', fv.faces');
    
    datalog.volume_fraction = sum(spheres_struct.vol)/totalVolume;
    
end

function totalVolume = renderSTL(fv)

%     figure(1)
    hold on;
     
    patch(fv,'FaceColor', [1 0 0], ...
         'EdgeColor', 'none', ...
         'FaceLighting', 'gouraud', ...
         'AmbientStrength', 0.85); 

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');
    alpha(0.3)                                                             
    axis('image');
%     view([-135 35]);
    view([180 0]);
    
    [totalVolume, totalArea] = stlVolume(fv.vertices', fv.faces');
    
end

function [data, fv_scaled] = STLPoints(data, fv_not_scaled, real_d_x_m)

    % Not scaled object
    tab = fv_not_scaled.vertices;
    
    x_all = tab(:,1);                       
    y_all = tab(:,2);
    z_all = tab(:,3);
    
    min_x = min(x_all);
    min_y = min(y_all);
    min_z = min(z_all);
    max_x = max(x_all);
    max_y = max(y_all);
    max_z = max(z_all);
    
    % Scaled image to real dimensions
    diff_pxl = abs(max_x-min_x);
%     c_factor = real_d_x_m/diff_pxl;
    c_factor = 1;
    
    fv_scaled.faces = fv_not_scaled.faces;
    fv_scaled.vertices = c_factor * fv_not_scaled.vertices;
    
    % Working with scaled object
    tab = [];
    x_all = [];
    y_all = [];
    z_all = [];
    min_x = [];
    max_x = [];
    min_y = [];
    max_y = [];
    min_z = [];
    max_z = [];
    
    tab = fv_scaled.vertices;
    
    % Total number of points in the STL file
    x_all = tab(:,1);                       
    y_all = tab(:,2);
    z_all = tab(:,3);
    
    % Number of points in the STL file
    N_x_tot = length(x_all);
    N_y_tot = N_x_tot;
    N_z_tot = N_x_tot;
    data.N_points_tot = N_x_tot;
    
    % Creation of a random permutation of the indexes of all the points
    N_rand_idxs_tot = randperm(N_x_tot);

    % Take the first N elements in the random vector
    N_rand_idxs = N_rand_idxs_tot(1:data.N_points_net);
    
    x_net = x_all(N_rand_idxs);
    y_net = y_all(N_rand_idxs);
    z_net = z_all(N_rand_idxs);
    
    hold off
    figure(3)
    plot3(x_net, y_net, z_net, '.')
    axis equal
    
    % Min points
    min_x = min(x_net);
    min_y = min(y_net);
    min_z = min(z_net);
    
    % Max points
    max_x = max(x_net);
    max_y = max(y_net);
    max_z = max(z_net);
    
    data.x_net = x_net;
    data.y_net = y_net;
    data.z_net = z_net;
    
    data.min_x = min_x;
    data.min_y = min_y;
    data.min_z = min_z;
    data.max_x = max_x;
    data.max_y = max_y;
    data.max_z = max_z;
    
end

function [x_c, y_c, z_c] = findCentroid(x, y, z)

    x_c = mean(x);
    y_c = mean(y);
    z_c = mean(z);

end

function [dist_net, pos_min, check, dists] = findClosestPoint(x, y, z, x_c, y_c, z_c, r_c)

    dists = sqrt( (x-x_c).^2 + (y-y_c).^2 + (z-z_c).^2);
    diff = dists-r_c;
    check = any((diff)<0);
    [dist_net, pos_min] = min(abs(diff));

end

function [points_x_c, points_y_c, points_z_c] = createSingleSphere(r, x_c, y_c, z_c)

    % General sphere generator for coating
    [x, y, z] = sphere(30);

    points_x_c = x * r + x_c;
    points_y_c = y * r + y_c;
    points_z_c = z * r + z_c;
    
end

function angle_deg = rad2deg(angle_rad)

    angle_deg = 180 .* angle_rad /pi;

end

function [totalVolume, totalArea] = stlVolume(p,t)

    % Given a surface triangulation, compute the volume enclosed using
    % divergence theorem.
    % Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
    % Input: p: (3xnPoints), t: (3xnTriangles)
    % Output: total volume enclosed, and total area of surface  
    % Author: K. Suresh; suresh@engr.wisc.edu

    % Compute the vectors d13 and d12
    d13 = [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
    d12 = [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
    cr = cross(d13,d12,1);                                                 % Cross-product (vectorized)
    area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);                     % Area of each triangle
    totalArea = sum(area);
    crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
    zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
    nz = -cr(3,:)./crNorm;                                                 % z component of normal for each triangle
    volume = area.*zMean.*nz;                                              % Contribution of each triangle
    totalVolume = sum(volume);                                             % Divergence theorem
    
end

function spheres_struct_light = filterFinalObject(spheres_struct, fv, N_spheres, alpha_transparency)

    x_coord = spheres_struct.x_c;
    y_coord = spheres_struct.y_c;
    z_coord = spheres_struct.z_c;
    P(:,1) = x_coord;
    P(:,2) = y_coord;
    P(:,3) = z_coord;
    k = boundary(P);
    all_indexes = [1:1:N_spheres];
    check_x = ismember(all_indexes,k(:,1));
    check_y = ismember(all_indexes,k(:,2));
    check_z = ismember(all_indexes,k(:,3));
    
    check_all = check_x+check_y+check_z;
    idxs_good = find(check_all>0);

    figure(5)
    
    for i=1:length(idxs_good)
        
        spheres_struct_light.r(i) = spheres_struct.r(idxs_good(i));
        spheres_struct_light.vol(i) = spheres_struct.vol(idxs_good(i));
        spheres_struct_light.x_c(i) = spheres_struct.x_c(idxs_good(i));
        spheres_struct_light.y_c(i) = spheres_struct.y_c(idxs_good(i));
        spheres_struct_light.z_c(i) = spheres_struct.z_c(idxs_good(i));
        spheres_struct_light.id(i) = i;
        
        [points_x_c, points_y_c, points_z_c] = createSingleSphere(spheres_struct_light.r(i), ...
                                                                  spheres_struct_light.x_c(i), ...
                                                                  spheres_struct_light.y_c(i), ...
                                                                  spheres_struct_light.z_c(i));
                    
         hl2 = surf(points_x_c, points_y_c, points_z_c);
         set(hl2,'FaceColor',[1 0 0], ...
                'FaceAlpha',alpha_transparency,'FaceLighting','gouraud','EdgeColor','none')
         alpha(hl2, alpha_transparency)
         hold on
         drawnow
         
    end
    
    axis equal
    renderSTL(fv);
    view([0 90])
    
   
end