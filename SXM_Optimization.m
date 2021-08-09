%% SXM Optimization

load('dp718_total_conf_refine_graindata.mat')
load('dp718_total_conf_refine_grainids_ele_new.mat')
load('dp718_total_conf_refine_grainids_ele.mat')
load('dp718_total_conf_refine_solution.mat')


%% Optimization of SXM c0_1 MM

% initialize constants
macro_stress = [0; 400; 0; 0; 0; 0]; % MPa
sx_mod_init = [259.0e3 176.0e3 106.0e3]; % MPa
nominal_strain_multiplier = 2.0;

max_iter = 10;
max_func_eval = 50;
options = optimset('MaxIter', max_iter, 'MaxFunEvals', max_func_eval, 'Display', 'iter');
options = optimset('MaxIter', max_iter, 'Display', 'iter');

[new_sx_mod_c0_1_MM, delta_stress] = fminsearch(@(SX_Mod)SXM_OptimizationMM_Strain_Avg(SX_Mod, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier), sx_mod_init, options);

% [vol_weighted_ave_microstructure_sigma] = SXM_Average_Stress(new_sx_mod_c0_1, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier);

new_sx_mod_c0_1_MM = new_sx_mod_c0_1_MM/1000.0;

fprintf('Optimized SXM MM c0_1 (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',new_sx_mod_c0_1_MM)


% fprintf('Volume Weighted Average Stress for Microstructure (MPa) \n')
% fprintf(' e_xx      e_yy      e_zz      e_yz      e_xz      e_xy \n')
% fprintf('%3.3f    %3.3f    %3.3f    %3.3f    %3.3f    %3.3f \n',vol_weighted_ave_microstructure_sigma)


%% Optimization of SXM c0_2 MM

% initialize constants
macro_stress = [0; 700; 0; 0; 0; 0]; % MPa
sx_mod_init = [259.0e3 176.0e3 106.0e3]; % MPa
nominal_strain_multiplier = 3.5;

max_iter = 50;
max_func_eval = 50;
options = optimset('MaxIter', max_iter, 'MaxFunEvals', max_func_eval, 'Display', 'iter');
options = optimset('MaxIter', max_iter, 'Display', 'iter');

[new_sx_mod_c0_2_MM, delta_stress] = fminsearch(@(SX_Mod)SXM_OptimizationMM_Strain_Avg(SX_Mod, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier), sx_mod_init, options);

% [vol_weighted_ave_microstructure_sigma] = SXM_Average_Stress(new_sx_mod_c0_2, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier);

new_sx_mod_c0_2_MM = new_sx_mod_c0_2_MM/1000.0;

fprintf('Optimized SXM MM c0_2 (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',new_sx_mod_c0_2_MM)


% fprintf('Volume Weighted Average Stress for Microstructure (MPa) \n')
% fprintf(' e_xx      e_yy      e_zz      e_yz      e_xz      e_xy \n')
% fprintf('%3.3f    %3.3f    %3.3f    %3.3f    %3.3f    %3.3f \n',vol_weighted_ave_microstructure_sigma)


%% SXM Optimization from FF Data

load('dp718_ff_rotation_mats.mat')
load('Sarah_ff_Values_corrected_conf_c0_1.mat')
load('dp718_total_conf_graindata.mat') % grain_volumes
ff = [ff_e_xx ff_e_yy ff_e_zz ff_e_yz ff_e_xz ff_e_xy]';
lit_val = [259.6 179.0 109.6]; % Literature values for c11, c12, c44


%% Optimization of SXM c0_1 FF

% initialize constants
macro_stress = [0; 400; 0; 0; 0; 0]; % MPa
sx_mod_init = [259.0e3 176.0e3 106.0e3]; % MPa
nominal_strain_multiplier = 1.0;

max_iter = 50;
max_func_eval = 50;
options = optimset('MaxIter', max_iter, 'MaxFunEvals', max_func_eval, 'Display', 'iter');
options = optimset('MaxIter', max_iter, 'Display', 'iter');

[new_sx_mod_c0_1_FF, delta_stress] = fminsearch(@(SX_Mod)SXM_OptimizationFF_Strain_Avg(SX_Mod, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier), sx_mod_init, options);

% [vol_weighted_ave_microstructure_sigma] = SXM_Average_Stress(new_sx_mod_c0_1_FF, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier);

new_sx_mod_c0_1_FF = new_sx_mod_c0_1_FF/1000.0;

fprintf('Optimized SXM FF c0_1 (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',new_sx_mod_c0_1_FF)


% fprintf('Volume Weighted Average Stress for Microstructure (MPa) \n')
% fprintf(' e_xx      e_yy      e_zz      e_yz      e_xz      e_xy \n')
% fprintf('%3.3f    %3.3f    %3.3f    %3.3f    %3.3f    %3.3f \n',vol_weighted_ave_microstructure_sigma)


%% Load c0_2 Values

load('Sarah_ff_Values_corrected_conf_c0_2.mat')
ff = [ff_e_xx ff_e_yy ff_e_zz ff_e_yz ff_e_xz ff_e_xy]';


%% Optimization of SXM c0_2 FF

% initialize constants
macro_stress = [0; 700; 0; 0; 0; 0]; % MPa
sx_mod_init = [259.0e3 176.0e3 106.0e3]; % MPa
nominal_strain_multiplier = 1.0;

max_iter = 50;
max_func_eval = 50;
options = optimset('MaxIter', max_iter, 'MaxFunEvals', max_func_eval, 'Display', 'iter');
options = optimset('MaxIter', max_iter, 'Display', 'iter');

[new_sx_mod_c0_2_FF, delta_stress] = fminsearch(@(SX_Mod)SXM_OptimizationFF_Strain_Avg(SX_Mod, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier), sx_mod_init, options);

% [vol_weighted_ave_microstructure_sigma] = SXM_Average_Stress(new_sx_mod_c0_2_FF, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier);

new_sx_mod_c0_2_FF = new_sx_mod_c0_2_FF/1000.0;

fprintf('Optimized SXM FF c0_2 (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',new_sx_mod_c0_2_FF)


% fprintf('Volume Weighted Average Stress for Microstructure (MPa) \n')
% fprintf(' e_xx      e_yy      e_zz      e_yz      e_xz      e_xy \n')
% fprintf('%3.3f    %3.3f    %3.3f    %3.3f    %3.3f    %3.3f \n',vol_weighted_ave_microstructure_sigma)


%% Average SXM (c11, c12, c44)

SXM_MM = [new_sx_mod_c0_1_MM; new_sx_mod_c0_2_MM];
ave_SXM_MM = mean(SXM_MM);

SXM_FF = [new_sx_mod_c0_1_FF; new_sx_mod_c0_2_FF];
ave_SXM_FF = mean(SXM_FF);

fprintf('Average SXM MM (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',ave_SXM_MM)

fprintf('Average SXM FF (GPa) \n')
fprintf(' c_11      c_12     c_44\n')
fprintf('%3.3f    %3.3f    %3.3f\n\n',ave_SXM_FF)

%% Plots

IDplot = (1:3); % c11, c12, c44
plot(IDplot, lit_val, 'k*', 'MarkerSize', 10)
hold on
grid on
plot(IDplot, ave_SXM_MM, 'm*', 'MarkerSize', 10)
plot(IDplot, ave_SXM_FF, 'g*', 'MarkerSize', 10)
xticklabels({'C11', ' ', 'C12', ' ', 'C44'})
title('Average Optimized Single Crystal Moduli') 
ylabel('Single Crystal Modulus (GPa)')
legend('Literature Values', 'MechMet', 'Far Field')








%% Optimization Function MM

function [optimize_val] = SXM_OptimizationMM(SX_Mod, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier)
% Number of grains
x = max(grains);

% Defining single crystal moduli and single crystal stiffness matrix in
% crystal coordinate system
c11 = SX_Mod(1); % MPa
c12 = SX_Mod(2); % MPa
c44 = SX_Mod(3); % MPa
sx_crys = [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];

          
% Gets sigma average (e_xx, e_yy, etc) for each grain using rotated
% stiffness matrix
ave_microstructure_sigma = zeros(x, 6);
          
for grain_ID = 1:x % rows
    rot = rotations(:, :, grain_ID); % finds rotation matrix of current grain ID
    sx = RotateStiffnessMatrix(sx_crys, rot); % creates stiffness matrix of current grain ID
    
    % calculates stresses of grain using rotated stiffness matrix
    epsilon_samp = epsilon_ave(:, grains == grain_ID) * nominal_strain_multiplier;
    sigma_samp = sx * epsilon_samp; % MPa
    ave_grain_sigma = mean(sigma_samp'); % average stress values for grain (6 strain tensor components)
    ave_microstructure_sigma(grain_ID, 1:6) = ave_grain_sigma;
end


% Gets volume weighted average of stresses (e_xx, e_yy, etc) for entire microstructure
vol_weighted_ave_microstructure_sigma = zeros(1, 6);
ID = (1:x);

for y = 1:6 % columns
    vol_sigma_sum = 0;
    for n = 1:length(ID) % rows equal to amount of grains
        if isnan(ave_microstructure_sigma(n, y))
            continue
        end
        vol_X_microstructure_sigma = grain_volumes(n) * ave_microstructure_sigma(n, y); 
        vol_sigma_sum = vol_sigma_sum + vol_X_microstructure_sigma; 
    end
    vol_weighted_ave_microstructure_sigma(1, y) = vol_sigma_sum/sum(grain_volumes); % MPa
end

temp = vol_weighted_ave_microstructure_sigma-macro_stress';
optimize_val = norm(vol_weighted_ave_microstructure_sigma-macro_stress'); % should be as close to 0 as possible

end

%% Optimization Function MM (Strain Averaging)

function [optimize_val] = SXM_OptimizationMM_Strain_Avg(SX_Mod, rotations, grains, epsilon_ave, grain_volumes, macro_stress, nominal_strain_multiplier)
% Number of grains
x = max(grains);

% Defining single crystal moduli and single crystal stiffness matrix in
% crystal coordinate system
c11 = SX_Mod(1); % MPa
c12 = SX_Mod(2); % MPa
c44 = SX_Mod(3); % MPa
sx_crys = [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];

          
% Gets sigma average (e_xx, e_yy, etc) for each grain using rotated
% stiffness matrix
ave_microstructure_sigma = zeros(x, 6);
          
for grain_ID = 1:x % rows
    rot = rotations(:, :, grain_ID); % finds rotation matrix of current grain ID
    sx = RotateStiffnessMatrix(sx_crys, rot); % creates stiffness matrix of current grain ID
    
    % calculates stresses of grain using rotated stiffness matrix
    epsilon_samp = epsilon_ave(:, grains == grain_ID) * nominal_strain_multiplier;
    ave_epsilon_samp = mean(epsilon_samp, 2); % 6 x 1 of avg epsilon values
    sigma_samp = sx * ave_epsilon_samp; % MPa 
    ave_grain_sigma = mean(sigma_samp, 2); % average stress values for grain (6 strain tensor components)
    ave_microstructure_sigma(grain_ID, 1:6) = ave_grain_sigma;
end


% Gets volume weighted average of stresses (e_xx, e_yy, etc) for entire microstructure
vol_weighted_ave_microstructure_sigma = zeros(1, 6);
ID = (1:x);

for y = 1:6 % columns
    vol_sigma_sum = 0;
    for n = 1:length(ID) % rows equal to amount of grains
        if isnan(ave_microstructure_sigma(n, y))
            continue
        end
        vol_X_microstructure_sigma = grain_volumes(n) * ave_microstructure_sigma(n, y); 
        vol_sigma_sum = vol_sigma_sum + vol_X_microstructure_sigma; 
    end
    vol_weighted_ave_microstructure_sigma(1, y) = vol_sigma_sum/sum(grain_volumes); % MPa
end

temp = vol_weighted_ave_microstructure_sigma-macro_stress';
optimize_val = norm(vol_weighted_ave_microstructure_sigma-macro_stress'); % should be as close to 0 as possible

end


%% Optimization Function FF

function [optimize_val] = SXM_OptimizationFF(SX_Mod, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier)
% Number of grains
x = max(ID);

% Defining single crystal moduli and single crystal stiffness matrix in
% crystal coordinate system
c11 = SX_Mod(1); % MPa
c12 = SX_Mod(2); % MPa
c44 = SX_Mod(3); % MPa
sx_crys = [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];

          
% Gets sigma average (e_xx, e_yy, etc) for each grain using rotated
% stiffness matrix
ave_microstructure_sigma = zeros(x, 6);
          
for grain_ID = 1:x % rows
    rot = ff_rot_mat(:, :, grain_ID); % finds rotation matrix of current grain ID
    sx = RotateStiffnessMatrix(sx_crys, rot); % creates stiffness matrix of current grain ID
    
    % calculates stresses of grain using rotated stiffness matrix
    epsilon_samp = ff(:, ID == grain_ID) * nominal_strain_multiplier;
    sigma_samp = sx * epsilon_samp; % MPa
    % ave_grain_sigma = mean(sigma_samp'); % average stress values for grain (6 strain tensor components)
    ave_microstructure_sigma(grain_ID, 1:6) = sigma_samp;
end


% Gets volume weighted average of stresses (e_xx, e_yy, etc) for entire microstructure
vol_weighted_ave_microstructure_sigma = zeros(1, 6);
ID = (1:x);

for y = 1:6 % columns
    vol_sigma_sum = 0;
    for n = 1:length(ID) % rows equal to amount of grains
        if isnan(ave_microstructure_sigma(n, y))
            continue
        end
        vol_X_microstructure_sigma = grain_volumes(n) * ave_microstructure_sigma(n, y); 
        vol_sigma_sum = vol_sigma_sum + vol_X_microstructure_sigma; 
    end
    vol_weighted_ave_microstructure_sigma(1, y) = vol_sigma_sum/sum(grain_volumes); % MPa
end

optimize_val = norm(vol_weighted_ave_microstructure_sigma-macro_stress'); % should be as close to 0 as possible

end


%% Optimization Function FF (Strain Averaging)

function [optimize_val] = SXM_OptimizationFF_Strain_Avg(SX_Mod, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier)
% Number of grains
x = max(ID);

% Defining single crystal moduli and single crystal stiffness matrix in
% crystal coordinate system
c11 = SX_Mod(1); % MPa
c12 = SX_Mod(2); % MPa
c44 = SX_Mod(3); % MPa
sx_crys = [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];

          
% Gets sigma average (e_xx, e_yy, etc) for each grain using rotated
% stiffness matrix
ave_microstructure_sigma = zeros(x, 6);
          
for grain_ID = 1:x % rows
    rot = ff_rot_mat(:, :, grain_ID); % finds rotation matrix of current grain ID
    sx = RotateStiffnessMatrix(sx_crys, rot); % creates stiffness matrix of current grain ID
    
    % calculates stresses of grain using rotated stiffness matrix
    epsilon_samp = ff(:, ID == grain_ID) * nominal_strain_multiplier;
    ave_epsilon_samp = mean(epsilon_samp, 2); % 6 x 1 of avg epsilon values
    sigma_samp = sx * ave_epsilon_samp; % MPa 
    ave_grain_sigma = mean(sigma_samp, 2); % average stress values for grain (6 strain tensor components)
    ave_microstructure_sigma(grain_ID, 1:6) = ave_grain_sigma;
end


% Gets volume weighted average of stresses (e_xx, e_yy, etc) for entire microstructure
vol_weighted_ave_microstructure_sigma = zeros(1, 6);
ID = (1:x);

for y = 1:6 % columns
    vol_sigma_sum = 0;
    for n = 1:length(ID) % rows equal to amount of grains
        if isnan(ave_microstructure_sigma(n, y))
            continue
        end
        vol_X_microstructure_sigma = grain_volumes(n) * ave_microstructure_sigma(n, y); 
        vol_sigma_sum = vol_sigma_sum + vol_X_microstructure_sigma; 
    end
    vol_weighted_ave_microstructure_sigma(1, y) = vol_sigma_sum/sum(grain_volumes); % MPa
end

optimize_val = norm(vol_weighted_ave_microstructure_sigma-macro_stress'); % should be as close to 0 as possible

end


%% SXM Average Stress Function FF

function [vol_weighted_ave_microstructure_sigma] = SXM_Average_Stress(SX_Mod, ff_rot_mat, ID, ff, grain_volumes, macro_stress, nominal_strain_multiplier)
% Number of grains
x = max(ID);

% Defining single crystal moduli and single crystal stiffness matrix in
% crystal coordinate system
c11 = SX_Mod(1); % MPa
c12 = SX_Mod(2); % MPa
c44 = SX_Mod(3); % MPa
sx_crys = [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];

          
% Gets sigma average (e_xx, e_yy, etc) for each grain using rotated
% stiffness matrix
ave_microstructure_sigma = zeros(x, 6);
          
for grain_ID = 1:x % rows
    rot = ff_rot_mat(:, :, grain_ID); % finds rotation matrix of current grain ID
    sx = RotateStiffnessMatrix(sx_crys, rot); % creates stiffness matrix of current grain ID
    
    % calculates stresses of grain using rotated stiffness matrix
    epsilon_samp = ff(:, ID == grain_ID) * nominal_strain_multiplier;
    sigma_samp = sx * epsilon_samp; % MPa
    % ave_grain_sigma = mean(sigma_samp'); % average stress values for grain (6 strain tensor components)
    ave_microstructure_sigma(grain_ID, 1:6) = sigma_samp;
end


% Gets volume weighted average of stresses (e_xx, e_yy, etc) for entire microstructure
vol_weighted_ave_microstructure_sigma = zeros(1, 6);
ID = (1:x);

for y = 1:6 % columns
    vol_sigma_sum = 0;
    for n = 1:length(ID) % rows equal to amount of grains
        if isnan(ave_microstructure_sigma(n, y))
            continue
        end
        vol_X_microstructure_sigma = grain_volumes(n) * ave_microstructure_sigma(n, y); 
        vol_sigma_sum = vol_sigma_sum + vol_X_microstructure_sigma; 
    end
    vol_weighted_ave_microstructure_sigma(1, y) = vol_sigma_sum/sum(grain_volumes); % MPa
end

end