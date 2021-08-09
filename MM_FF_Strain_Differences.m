%% Find stress/strain from MechMet, difference between MechMet and FF, Weighted Averages 
%% c0_1 400 MPa

clc, clear all

load('dp718_total_conf_refine_grainids_ele.mat')
load('dp718_total_conf_refine_grainids_ele_new.mat')
load('dp718_total_conf_refine_solution.mat')
load('dp718_total_conf_refine_graindata.mat')
load('Sarah_ff_Values_corrected_conf_c0_1.mat')

ff=[ff_e_xx ff_e_yy ff_e_zz ff_e_yz ff_e_xz ff_e_xy];

% Convert MechMet values to appropriate strain
% c0_1 = 1.994
s = 1.994;
simga_ave = sigma_ave*s;
epsilon_ave = epsilon_ave*s;
simga = sigma*s;
epsilon = epsilon*s;

% Find avg of stress/strain for each grain from MechMet

x=max(grains); 

grain_avg_stress = zeros(x,6); 
grain_avg_strain = zeros(x,6); 

grain_stress = zeros(x,6);
grain_strain = zeros(x,6);

for grainid = 1:x   % rows
    grain_index = find(grain4np == grainid);
    grain_index_ave = find(grains == grainid);
    
    for n = 1:6 % e_xx, e_yy, e_zz, e_yz, e_xz, e_xy (columns)
        stress_sig_ave = mean(sigma_ave(n,grain_index_ave));
        strain_eps_ave = mean(epsilon_ave(n,grain_index_ave));
        
        stress_ave = mean(epsilon(n,grain_index));
        strain_ave = mean(sigma(n,grain_index));
        
        grain_avg_stress(grainid,n) = stress_sig_ave;
        grain_avg_strain(grainid,n) = strain_eps_ave;
        
        grain_stress(grainid,n) = stress_ave;
        grain_strain(grainid,n) = strain_ave;
    end
end

% MechMet volume weighted average of average strain

volume_weightedMM=zeros(1,6);

for y= 1:6
    volume_strain_sumMM=0;
    for n= 1:length(ID)
        if isnan(grain_avg_strain(n,y))
            continue
        end
        volumeXstrainMM = grain_volumes(n)*grain_avg_strain(n,y);
        volume_strain_sumMM = volume_strain_sumMM+volumeXstrainMM;
    end
    volume_weightedMM(1,y) = volume_strain_sumMM/sum(grain_volumes);
end

% FF volume weighted average of average strain

volume_weightedFF=zeros(1,6);

for z= 1:6
    volume_strain_sumFF=0;
    for n= 1:length(ID)
        if isnan(ff(n,z))
            continue
        end
        volumeXstrainFF = grain_volumes(n)*ff(n,z);
        volume_strain_sumFF = volume_strain_sumFF+volumeXstrainFF;
    end
    volume_weightedFF(1,z) = volume_strain_sumFF/sum(grain_volumes);
end

% Difference between ff and MechMet

ff_threshold = ff(((grain_volumes >= 2e-4)' & (completeness >= 0.85) & (chi2 <= 1e-2)), :);
grain_avg_strain_threshold = grain_avg_strain(((grain_volumes >= 2e-4)' & (completeness >= 0.85) & (chi2 <= 1e-2)), :);

diff_exx=ff_threshold(:,1)-grain_avg_strain_threshold(:,1);
diff_eyy=ff_threshold(:,2)-grain_avg_strain_threshold(:,2);
diff_ezz=ff_threshold(:,3)-grain_avg_strain_threshold(:,3);
diff_eyz=ff_threshold(:,4)-grain_avg_strain_threshold(:,4);
diff_exz=ff_threshold(:,5)-grain_avg_strain_threshold(:,5);
diff_exy=ff_threshold(:,6)-grain_avg_strain_threshold(:,6);

% Entire microstructure average

avg_mechmet=mean(grain_avg_strain,'omitnan');
avg_ff=mean(ff,'omitnan');

ffmax=max(ff);
ffmin=min(ff);

diff_ff_mechmet=[diff_exx diff_eyy diff_ezz diff_eyz diff_exz diff_exy];

avg_diff=mean(diff_ff_mechmet,'omitnan');

IDavg=(1:6);

figure(1)
plot(IDavg, volume_weightedMM, 'g*', 'MarkerSize', 10)
hold on
plot(IDavg, volume_weightedFF, 'b*', 'MarkerSize', 10)
grid on 
title('Average Strain Data at 400 MPa')
xticklabels({'e xx','e yy','e zz','e yz','e xz','e xy'})
ylabel('Average Strain')
legend('MechMet Volume-Weighted', 'Far Field Volume-Weighted')

% Standard Deviation

std_diff_exx=std(diff_exx,'omitnan');
std_diff_eyy=std(diff_eyy,'omitnan');
std_diff_ezz=std(diff_ezz,'omitnan');
std_diff_eyz=std(diff_eyz,'omitnan');
std_diff_exz=std(diff_exz,'omitnan');
std_diff_exy=std(diff_exy,'omitnan');
std=[std_diff_exx std_diff_eyy std_diff_ezz std_diff_eyz std_diff_exz std_diff_exy];


%%
fprintf('Microstructure Average c0 1, 400 MPa \n')
fprintf('Far Field \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_ff)
fprintf('Far Field Weighted \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',volume_weightedFF)
fprintf(' \n')
fprintf('MechMet \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_mechmet)
fprintf('MechMet Weighted \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',volume_weightedMM)

fprintf(' \n')
fprintf('Average Difference between MechMet and FF \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_diff)

fprintf(' \n')
fprintf('Far Field max \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',ffmax)
fprintf('Far Field min \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',ffmin)

fprintf(' \n')
fprintf('Standard Deviation of Difference \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n \n \n',std)


%%  Histogram with Threshold c0_1

threshold=2e-4;

figure(2)
sgtitle('Difference between MechMet and FF Strain Values at 400 MPa')
subplot(2,3,1)
h1 = histogram(diff_exx);
h1.BinWidth = 1e-4;
title('e xx')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,2)
h2 = histogram(diff_eyy);
h2.BinWidth = 1e-4;
title('e yy')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,3)
h3 = histogram(diff_ezz);
h3.BinWidth = 1e-4;
title('e zz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,4)
h4 = histogram(diff_eyz);
h4.BinWidth = 1e-4;
title('e yz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,5)
h5 = histogram(diff_exz);
h5.BinWidth = 1e-4;
title('e xz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,6)
h6 = histogram(diff_exy);
h6.BinWidth = 1e-4;
title('e xy')
xlim([-3e-3 3e-3])
ylim([0 110])


%% c0_2 400 MPa

clear all

load('dp718_total_conf_refine_grainids_ele.mat')
load('dp718_total_conf_refine_grainids_ele_new.mat')
load('dp718_total_conf_refine_solution.mat')
load('dp718_total_conf_refine_graindata.mat')
load('Sarah_ff_Values_corrected_conf_c0_2.mat')

ff=[ff_e_xx ff_e_yy ff_e_zz ff_e_yz ff_e_xz ff_e_xy];

% Convert MechMet values to appropriate strain
% c0_2 = 3.482
s = 3.482;
simga_ave = sigma_ave*s;
epsilon_ave = epsilon_ave*s;
simga = sigma*s;
epsilon = epsilon*s;

% Find avg of stress/strain for each grain from MechMet

x=max(grains); 

grain_avg_stress = zeros(x,6); 
grain_avg_strain = zeros(x,6); 

grain_stress = zeros(x,6);
grain_strain = zeros(x,6);

for grainid = 1:x   % rows
    grain_index = find(grain4np == grainid);
    grain_index_ave = find(grains == grainid);
    
    for n = 1:6 % e_xx, e_yy, e_zz, e_yz, e_xz, e_xy (columns)
        stress_sig_ave = mean(sigma_ave(n,grain_index_ave));
        strain_eps_ave = mean(epsilon_ave(n,grain_index_ave));
        
        stress_ave = mean(epsilon(n,grain_index));
        strain_ave = mean(sigma(n,grain_index));
        
        grain_avg_stress(grainid,n) = stress_sig_ave;
        grain_avg_strain(grainid,n) = strain_eps_ave;
        
        grain_stress(grainid,n) = stress_ave;
        grain_strain(grainid,n) = strain_ave;
    end
end

% MechMet volume weighted average of average strain

volume_weightedMM=zeros(1,6);

for y= 1:6
    volume_strain_sumMM=0;
    for n= 1:length(ID)
        if isnan(grain_avg_strain(n,y))
            continue
        end
        volumeXstrainMM = grain_volumes(n)*grain_avg_strain(n,y);
        volume_strain_sumMM = volume_strain_sumMM+volumeXstrainMM;
    end
    volume_weightedMM(1,y) = volume_strain_sumMM/sum(grain_volumes);
end

% FF volume weighted average of average strain

volume_weightedFF=zeros(1,6);

for z= 1:6
    volume_strain_sumFF=0;
    for n= 1:length(ID)
        if isnan(ff(n,z))
            continue
        end
        volumeXstrainFF = grain_volumes(n)*ff(n,z);
        volume_strain_sumFF = volume_strain_sumFF+volumeXstrainFF;
    end
    volume_weightedFF(1,z) = volume_strain_sumFF/sum(grain_volumes);
end

% Difference between ff and MechMet

ff_threshold = ff(((grain_volumes >= 2e-4)' & (completeness >= 0.85) & (chi2 <= 1e-2)), :);
grain_avg_strain_threshold = grain_avg_strain(((grain_volumes >= 2e-4)' & (completeness >= 0.85) & (chi2 <= 1e-2)), :);

diff_exx=ff_threshold(:,1)-grain_avg_strain_threshold(:,1);
diff_eyy=ff_threshold(:,2)-grain_avg_strain_threshold(:,2);
diff_ezz=ff_threshold(:,3)-grain_avg_strain_threshold(:,3);
diff_eyz=ff_threshold(:,4)-grain_avg_strain_threshold(:,4);
diff_exz=ff_threshold(:,5)-grain_avg_strain_threshold(:,5);
diff_exy=ff_threshold(:,6)-grain_avg_strain_threshold(:,6);

% Entire microstructure average

avg_mechmet=mean(grain_avg_strain,'omitnan');
avg_ff=mean(ff,'omitnan');

ffmax=max(ff);
ffmin=min(ff);

diff_ff_mechmet=[diff_exx diff_eyy diff_ezz diff_eyz diff_exz diff_exy];

avg_diff=mean(diff_ff_mechmet,'omitnan');

IDavg=(1:6);

figure(3)
plot(IDavg, volume_weightedMM, 'g*', 'MarkerSize', 10)
hold on
plot(IDavg, volume_weightedFF, 'b*', 'MarkerSize', 10)
grid on 
title('Average Strain Data at 700 MPa')
xticklabels({'e xx','e yy','e zz','e yz','e xz','e xy'})
ylabel('Average Strain')
legend('MechMet Volume-Weighted', 'Far Field Volume-Weighted')

% Standard Deviation

std_diff_exx=std(diff_exx,'omitnan');
std_diff_eyy=std(diff_eyy,'omitnan');
std_diff_ezz=std(diff_ezz,'omitnan');
std_diff_eyz=std(diff_eyz,'omitnan');
std_diff_exz=std(diff_exz,'omitnan');
std_diff_exy=std(diff_exy,'omitnan');
std=[std_diff_exx std_diff_eyy std_diff_ezz std_diff_eyz std_diff_exz std_diff_exy];

%%
fprintf('Microstructure Average c0 2, 700 MPa \n')
fprintf('Far Field \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_ff)
fprintf('Far Field Weighted \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',volume_weightedFF)
fprintf(' \n')
fprintf('MechMet \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_mechmet)
fprintf('MechMet Weighted \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',volume_weightedMM)

fprintf(' \n')
fprintf('Average Difference between MechMet and FF \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',avg_diff)

fprintf(' \n')
fprintf('Far Field max \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',ffmax)
fprintf('Far Field min \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',ffmin)

fprintf(' \n')
fprintf('Standard Deviation of Difference \n')
fprintf('   e xx       e yy       e zz       e yz       e xz       e xy\n')
fprintf('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n',std)


%%  Histogram with Threshold c0_2

threshold=2e-4;

figure(4)
sgtitle('Difference between MechMet and FF Strain Values at 700 MPa')
subplot(2,3,1)
h1 = histogram(diff_exx);
h1.BinWidth = 1e-4;
title('e xx')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,2)
h2 = histogram(diff_eyy);
h2.BinWidth = 1e-4;
title('e yy')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,3)
h3 = histogram(diff_ezz);
h3.BinWidth = 1e-4;
title('e zz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,4)
h4 = histogram(diff_eyz);
h4.BinWidth = 1e-4;
title('e yz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,5)
h5 = histogram(diff_exz);
h5.BinWidth = 1e-4;
title('e xz')
xlim([-3e-3 3e-3])
ylim([0 110])

subplot(2,3,6)
h6 = histogram(diff_exy);
h6.BinWidth = 1e-4;
title('e xy')
xlim([-3e-3 3e-3])
ylim([0 110])

%% Histogram of Grain Size Distribution

% figure(3)
% histogram(grain_volumes,250)
% title('Grain Size Distribution of Microstructure')
% xlabel('Grain Volume')
% ylabel('Frequency')








