
% Important: 
% 1) copy the present script inside the folder 'SCARLET_v1/source_code'
% 2) copy the file "SCARLET_v1_GMDplotFig18.m" inside the folder
% "source_code"
% 3) set the Matlab path as 'SCARLET_v1/source_code'
% 4) Uncomment the part of the text that generates the figure you are
% interested in (uncomment between signs "...%%%%%%%%%%%%..."). A part is
% dedicated to reproduce ellipsoids; the other one to reproduce spheres.
% 5) Set "sort_boole = 0;" for sequential displacement; "sort_boole = 1;"
% for sequential one
% 6) type the name of the script in the Matlab Command Window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of plots reported in Fig. 18

close all;
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AGG EYJA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sort_boole = 0;

% -------------------------- %
% Import the STL shapes already filled with spheres: the structure of the
% input file for this test is made by the coarse ellipsoid at the center
% acting as a core; then the coarse ellipsoid and the fine one in position
% 1 and 2 respectively.

% PLEASE NOTE: it is the same file of Fig.16

load('../STL_examples/Fig15/aggr_27_sk21.mat');                            % Now we have agg_27_sk21 in the workspace

% Core size
core_size = 40; 

% Sizes of the coating for aggregate EJ15a from Bonaddona et al 2011 (in mm)
vector_sizes = [0.02096272 0.031734272 0.0092243309 0.0078826909 0.0075748663 ...
0.0145657394 0.0128288191 0.0139573372 0.0169791742 0.0187122767 0.0188315598 ...
0.0127302589 0.0225515009 0.0171921401 0.0105545474 0.0185997534 0.0284263671 ...
0.0316459494 0.0144014773 0.0157742554 0.0204969067 0.0258207491 0.0201665107 ...
0.0147566502 0.0182810089 0.0231395792 0.015296817 0.0111100803 0.007775439 ...
0.01036725 0.0177528108 0.0143137545 0.0138364929 0.017094179 0.0194385951 ...
0.0085307014 0.0123847811 0.0100380486 0.0256957877 0.0171677004 0.0140173687 ...
0.0190898913 0.0071958764 0.009869331 0.011982753 0.0167049632 0.0172490266 ...
0.0202427025 0.0134778422 0.0117706533 0.0128941057 0.0284312894 0.0188092503 ...
0.0086123375 0.0124298951 0.0193085671 0.0107124786 0.0177764462 0.0180963451 ...
0.0089470701 0.0192868128 0.0192068293 0.012686208 0.0296835153 0.0214901755 ...
0.0231999808 0.0125976492 0.0126641276 0.0091634416 0.0266528568 0.0239245976 ...
0.0126420086 0.0166041281 0.0122141038 0.0106074484 0.0147091548 0.0127741578 ...
0.0089783036 0.01711872 0.0186523474 0.0124073586 0.0167300789 0.014935752 ...
0.0212017035 0.0112353401 0.0090559056 0.016813522 0.0288418 0.0165788254 ...
0.008931416 0.0125196457 0.0128288191 0.010093657 0.0110215539 0.0152418256 ...
0.0169791742 0.0147566502 0.009434339 0.012396075 0.0084482765 0.0202081052 ...
0.0197599157 0.0231214269 0.013205109 0.0208824615 0.015836236 0.0242614974 ...
0.0188909172 0.0103537433 0.0105677975 0.0160730696 0.0155328663 0.0131200532 ...
0.0124636286 0.0139473076 0.0191776608 0.0346397474 0.0143235275 0.014097016 ...
0.0093149151 0.0150012017 0.0112104024 0.0247526026 0.0085142817 0.022975703 ...
0.0144983285 0.0133840622 0.017305726 0.0242095259 0.0178000501 0.0175705638 ...
0.0081960298 0.0062149501 0.0121681845 0.0128615013 0.0059385923 0.0076849178 ...
0.0087892668 0.0219795272 0.0146615056 0.0200481903 0.0185997534 0.0092697302 ...
0.0056487189 0.0139071122 0.0117468521 0.0171595473 0.0059621097 0.0076300863 ...
0.0076484127 0.0137858277 0.0069384464 0.0074067193 0.0087413652 0.0090559056 ...
0.0147756039 0.0090249407 0.0197386555 0.0159419289 0.0099540473 0.0158715465 ...
0.0114817607 0.0272397331 0.0128832484 0.0069986896 0.0283721605 0.0241110541 ...
0.0113345534 0.0151589666 0.0110974778 0.0080929316 0.0193520091 0.0071568772 ...
0.0052905155 0.0114939419 0.0092091411 0.0080234608 0.0184562419 0.0113221951 ...
0.0143625575 0.0297117897 0.0086123375 0.0091328483 0.0165619316 0.0162290374 ...
0.0150849273 0.0055487418 0.0184030875 0.0072152969 0.0219412911 0.0071764033 ...
0.0095229257 0.0064795317 0.0288029568 0.0150663628 0.0155328663 0.0109961313 ...
0.0146424042 0.0147945376 0.0067751607 0.0149263798 0.0312274867 0.0154967841 ...
0.0219604189 0.0191776608 0.0078826909 0.0074255884];

vettore_taglie_2 = vector_sizes*1000;           
    
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
N_repetitions = 1;

for i=1:N_repetitions
    [output_st, vector_collisions] = SCARLET_v1_GMDplotFig18(agg_27_sk21, vettore_taglie, ...
                                                             core_size, sort_boole, ...
                                                             cone, rays, euler);
    matrix_data(i,:) = vector_collisions;
end
  

for i=1:length(vector_collisions)
    
    y(i) = mean(matrix_data(:,i));
    dy(i) = std(matrix_data(:,i));
    
end

x = [1:1:length(vector_collisions)];

x_EYJA = x; 
y_EYJA = y;
dy_EYJA = dy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AGG 27sk21 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
N_repetitions = 1;

for i=1:N_repetitions
    [output_st, vector_collisions] = SCARLET_v1_GMDplotFig18(agg_27_sk21, vettore_taglie, ...
                                                             core_size, sort_boole, ...
                                                             cone, rays, euler);
    matrix_data(i,:) = vector_collisions;
end
  

for i=1:length(vector_collisions)
    
    y(i) = mean(matrix_data(:,i));
    dy(i) = std(matrix_data(:,i));
    
end

x = [1:1:length(vector_collisions)];

x_Sakura = x; 
y_Sakura = y;
dy_Sakura = dy;

figure(1)
shadedErrorBar(x_EYJA, y_EYJA, dy_EYJA,'lineProps','r')
hold on
shadedErrorBar(x_Sakura, y_Sakura, dy_Sakura,'lineProps','b')
ylim([0 0.9]);


