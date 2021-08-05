clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
% -------------------------------------------------------------------------
% 1.- Define & Generate target trajetories
% -------------------------------------------------------------------------
dt = (1e-1);
%target = oval_trayectory();
%target = car_up_left();
target = car_up_left_v2();
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
bss = BSs();
bss = bss.gen_pilot_tx();
fc = FC();

% trackers:
lkf = LKF_xy_based();
ekf = EKF_range_based();
ukf = UKF_range_based();

% params and hist
scene = Params.get_scene();
N_t = size(target.t_vect,2);
deltas_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);

xy_ekf_hist = zeros(4,N_t);
eig_Pekf_est_hist = zeros(4,N_t);
eig_Pekf_pred_hist = zeros(4,N_t);

xy_ukf_hist = zeros(4,N_t);
eig_Pukf_est_hist = zeros(4,N_t);
eig_Pukf_pred_hist = zeros(4,N_t);

xy_lkf_hist = zeros(4,N_t);
eig_Plkf_est_hist = zeros(4,N_t);
eig_Plkf_pred_hist = zeros(4,N_t);

for t_idx = 1:N_t
    target.history(:,t_idx);
    
    % noise_type: SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx));
    [~, deltas_mean, deltas_var,act_bss] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
    
    if t_idx == 1
        lkf = lkf.set_x0(xy_toa);
        ekf = ekf.set_x0(xy_toa);
        ukf = ukf.set_x0(xy_toa);
    else
        lkf = lkf.predict(dt);
        lkf = lkf.correct(xy_toa, varxy_toa);
        
        ekf = ekf.predict(dt);
        ekf = ekf.correct(deltas_mean,deltas_var);
        
        ukf = ukf.predict(dt);
        ukf = ukf.correct(deltas_mean,deltas_var);
    end
    
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    xy_toa_hist(:,t_idx) = xy_toa;
    % ekf
    xy_ekf_hist(:,t_idx) = ekf.x_est;
    eig_Pekf_est_hist(:,t_idx) = ekf.eig_P_est;
    eig_Pekf_pred_hist(:,t_idx) = ekf.eig_P_pred;
    % ukf
    xy_ukf_hist(:,t_idx) = ukf.x_est;
    eig_Pukf_est_hist(:,t_idx) = ukf.eig_P_est;
    eig_Pukf_pred_hist(:,t_idx) = ukf.eig_P_pred;
    % lkf
    xy_lkf_hist(:,t_idx) = lkf.x_est;
    eig_Plkf_est_hist(:,t_idx) = lkf.eig_P_est;
    eig_Plkf_pred_hist(:,t_idx) = lkf.eig_P_pred;
end

fig1 = figure('Position',[1925 847 560 420]);
show_target_and_xy_toa_in_state_space(fig1,target,xy_toa_hist);

fig2 = figure('Position',[2577 858 560 420]);
show_target_and_tracker(fig2, target, xy_ekf_hist,eig_Pekf_est_hist,eig_Pekf_pred_hist, 'EKF');

fig3 = figure('Position',[3158 562 560 420]);
show_target_and_tracker(fig3, target, xy_ukf_hist,eig_Pukf_est_hist,eig_Pukf_pred_hist, 'UKF');

fig5 = figure('Position',[2554 52 560 420]);
show_target_and_tracker(fig5, target, xy_lkf_hist,eig_Plkf_est_hist,eig_Plkf_pred_hist, 'LKF');

fig4 = figure('Position',[1943 38 560 420]);
show_target_meas_xy_toa_in_meas_space(fig4,target,xy_toa_hist, deltas_hist);

fig6 = figure('Position',[3159 42 560 420]);
kfs = {xy_lkf_hist,xy_ekf_hist,xy_ukf_hist};
kfs_labels = {'LKF','EKF','UKF'};
show_rmse_comparison(fig6,target,kfs, kfs_labels);