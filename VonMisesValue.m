%% Calculates max Von Mises stresses and finds grains with highest stress (c0_2 700 MPa)

load('dp718_total_conf_refine_700_solution_vm.mat') % vm_result
load('dp718_total_conf_refine_700_grainids_ele_new.mat') % grain4np
load('dp718_total_conf_refine_700_graindata.mat') % grain_volumes

x = max(grain4np);
prob_grains = zeros(x, 1);
prob_grainids = zeros(x, 1);
prob_grains_count = 0;

% finds grain IDs of grains experiencing Von Mises stresses greater than 1100 MPa
for grainid = 1:x   % rows
    grain_index = find(grain4np == grainid);
    grain_vol = grain_volumes(grainid);
    vm_grain = vm_result(grain_index);
    greater_than1100 = vm_grain >= 1100;
    if (max(greater_than1100) == 1) & (grain_vol >= 2e-4);
        prob_grains_count = prob_grains_count + 1;
        prob_grainids(prob_grains_count) = grainid;
        prob_grains(grainid) = 1;
        % fprintf('%d \n', grainid);
    else
        prob_grains(grainid) = 0;
    end
end

% list of grains with volume >= 2e-4 and Von Mises stress >= 1100 MPa:
prob_grainids = prob_grainids(1: prob_grains_count);


% finds individual Von Mises stresses within a specified grain and max stress
% within that grain (grain 44 used as example)
grainid_single = 44;   % rows
grain_index_single = find(grain4np == grainid_single);
vm_grain_single = vm_result(grain_index_single);
max_single_grain = max(vm_grain_single);


%% Calculates max Von Mises stresses and finds grains with highest stress (c0_1 400 MPa)

load('dp718_total_conf_refine_400_solution_vm.mat') % vm_result
load('dp718_total_conf_refine_400_grainids_ele_new.mat') % grain4np
load('dp718_total_conf_refine_400_graindata.mat') % grain_volumes

x = max(grain4np);
prob_grains = zeros(x, 1);
prob_grainids = zeros(x, 1);
prob_grains_count = 0;

% finds grain IDs of grains experiencing Von Mises stresses greater than 1100 MPa
for grainid = 1:x   % rows
    grain_index = find(grain4np == grainid);
    grain_vol = grain_volumes(grainid);
    vm_grain = vm_result(grain_index);
    greater_than1100 = vm_grain >= 1100;
    if (max(greater_than1100) == 1) & (grain_vol >= 2e-4);
        prob_grains_count = prob_grains_count + 1;
        prob_grainids(prob_grains_count) = grainid;
        prob_grains(grainid) = 1;
        % fprintf('%d \n', grainid);
    else
        prob_grains(grainid) = 0;
    end
end

% list of grains with volume >= 2e-4 and Von Mises stress >= 1100 MPa:
prob_grainids = prob_grainids(1: prob_grains_count);


% finds individual Von Mises stresses within a specified grain and max stress
% within that grain (grain 44 used as example)
grainid_single = 44;   % rows
grain_index_single = find(grain4np == grainid_single);
vm_grain_single = vm_result(grain_index_single);
max_single_grain = max(vm_grain_single);













%%

% function [VonMisesValue] = VonMisesValue(sigma, grain4np)
% 
% % Calculates Von Mises stress given stress tensor values (sigma from
% % MechMet.m)
% 
% x= max(grain4np);
% sigma_ = sigma';
% 
% for grainid = 1:x   % rows
%     g_i = find(grain4np == grainid); % grain index
% 
%     s11 = sigma_(g_i, 1);
%     s22 = sigma_(g_i, 2);
%     s33 = sigma_(g_i, 3);
%     s23 = sigma_(g_i, 4);
%     s31 = sigma_(g_i, 5);
%     s12 = sigma_(g_i, 6);
%     
%     VM = sqrt(((s11-s22).^2+(s22-s33).^2+(s33-s11).^2+6.*(s12.^2+s23.^2+s31.^2))./2);
%     VM_avg = mean(VM);
%     
%     VonMisesValue(grainid) = VM_avg;
% end
% end
