clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
%rng(1)
% -------------------------------------------------------------------------
% 1.- Define & Generate target trajetories
% -------------------------------------------------------------------------
dt = (1e-1);
target = oval_trayectory();
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
bss = BSs();
bss = bss.gen_pilot_tx();
fc = FC();

% trackers:
iekf = IEKF_bsbnd_based();
dpekf = DPEKF_bsbnd_based();

% params and hist
scene = Params.get_scene();
N_t = size(target.t_vect,2);
deltas_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);

xy_iekf_hist = zeros(4,N_t);
eig_Piekf_est_hist = zeros(4,N_t);
eig_Piekf_pred_hist = zeros(4,N_t);

xy_dpekf_hist = zeros(4,N_t);
eig_Pdpekf_est_hist = zeros(2,N_t);
eig_Pdpekf_pred_hist = zeros(2,N_t);

for t_idx = 1:N_t
    target.history(:,t_idx)
    xy_true = [1 0 0 0; 0 0 1 0]*target.history(:,t_idx);
    % noise_type: SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx));
    [~, deltas_mean, deltas_var,act_bss] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
    [xy_bsbnd, varxy_bsbnd] = fc.dpe(bss.pilot_rx,bss.vars, bss.pilot_tx);
    
    if t_idx == 1
        iekf = iekf.set_x0(xy_true);
        dpekf = dpekf.set_x0(xy_true);
    else
        iekf = iekf.predict(dt);
        iekf = iekf.correct(squeeze(bss.pilot_rx),squeeze(bss.vars), squeeze(bss.pilot_tx));
        
        dpekf = dpekf.predict(dt);
        dpekf = dpekf.correct(bss,dt);
    end
    
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    xy_toa_hist(:,t_idx) = xy_bsbnd;
    % iekf
    xy_iekf_hist(:,t_idx) = iekf.x_est;
    eig_Piekf_est_hist(:,t_idx) = iekf.eig_P_est;
    eig_Piekf_pred_hist(:,t_idx) = iekf.eig_P_pred;
    % dpekf
    xy_dpekf_hist(:,t_idx) = [dpekf.x_est(1) dpekf.u(1) dpekf.x_est(2) dpekf.u(2)].';
    eig_Pdpekf_est_hist(:,t_idx) = dpekf.eig_P_est;
    eig_Pdpekf_pred_hist(:,t_idx) = dpekf.eig_P_pred;
end

fig1 = figure('Position',[1925 847 560 420]);
show_target_and_xy_toa_in_state_space(fig1,target,xy_toa_hist);

fig2 = figure('Position',[2577 858 560 420]);
show_target_and_tracker(fig2, target, xy_iekf_hist,eig_Piekf_est_hist,eig_Piekf_pred_hist, 'IEKF');

fig3 = figure('Position',[3211 571 570 413]);
show_target_and_tracker(fig3, target, xy_dpekf_hist,eig_Pdpekf_est_hist,eig_Pdpekf_pred_hist, 'DPEKF');

fig4 = figure('Position',[1943 38 560 420]);
show_target_meas_xy_toa_in_meas_space(fig4,target,xy_toa_hist, deltas_hist);

fig6 = figure('Position',[3159 42 560 420]);
kfs = {xy_dpekf_hist,xy_iekf_hist};
kfs_labels = {'DPEKF','IEKF'};
show_rmse_comparison(fig6,target,kfs, kfs_labels);