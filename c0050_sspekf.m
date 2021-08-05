clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
rng(2)
% -------------------------------------------------------------------------
% 1.- Define & Generate target trajetories
% -------------------------------------------------------------------------
dt = (1e-1);
%target = oval_trayectory();
target = car_up_left_v2();
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
bss = BSs();
bss = bss.gen_pilot_tx();
fc = FC();

% trackers:
dpekf = DPEKF_bsbnd_based();
%sspekf = SSPEKF_bsbnd_based_comp_sspe();
sspekf = SSPEKF_bsbnd_based_comp_nhop();
sspekf_lkf = SSPEKF_bsbnd_based_comp_lkf();

% params and hist
scene = Params.get_scene();
N_t = size(target.t_vect,2);
deltas_hist = zeros(scene.N_bs,N_t);

xy_bsbnd_hist = zeros(2,N_t);

xy_dpekf_hist = zeros(4,N_t);
eig_Pdpekf_est_hist = zeros(2,N_t);
eig_Pdpekf_pred_hist = zeros(2,N_t);


xy_sspekf_hist = zeros(4,scene.N_bs,N_t);
eig_Psspekf_est_hist = zeros(2,scene.N_bs,N_t);
eig_Psspekf_pred_hist = zeros(2,scene.N_bs,N_t);

for t_idx = 1:N_t
    target.history(:,t_idx)
    xy_true = [1 0 0 0; 0 0 1 0]*target.history(:,t_idx);
    
    % bsbnd signals + dpe(no prior)
    bss = bss.channel_propagation(target.history(:,t_idx));
    [xy_bsbnd, varxy_bsbnd] = fc.dpe(bss.pilot_rx, bss.vars, bss.pilot_tx);
    
    xy_0 = xy_true + 1* randn(2,1);
    xy_0_v = repmat(xy_true,1,scene.N_bs) + 1* randn(2,scene.N_bs);
    if t_idx == 1
        dpekf = dpekf.set_x0(xy_0);
        sspekf = sspekf.set_x0(xy_0_v);
    else
        dpekf = dpekf.predict(dt);
        dpekf = dpekf.correct(bss,dt);
        
        sspekf = sspekf.predict(dt);
        sspekf = sspekf.correct(bss,dt,t_idx);
    end
    
    % ______________________________________________
    % save for history:
    % bsbnd:
    xy_bsbnd_hist(:,t_idx) = xy_bsbnd;
    
    % dpekf
    xy_dpekf_hist(:,t_idx) = [dpekf.x_est(1) dpekf.u(1) dpekf.x_est(2) dpekf.u(2)].';
    eig_Pdpekf_est_hist(:,t_idx) = dpekf.eig_P_est;
    eig_Pdpekf_pred_hist(:,t_idx) = dpekf.eig_P_pred;
    
    % sspekf
    xy_sspekf_hist(:,:,t_idx) = sspekf.x_plot;
    eig_Psspekf_est_hist(:,:,t_idx) = sspekf.eig_P_est_v;
    eig_Psspekf_pred_hist(:,:,t_idx) = sspekf.eig_P_pred_v;
    
end

fig1 = figure('Position',[25   547   560   420]);
show_target_and_xy_toa_in_state_space(fig1,target,xy_bsbnd_hist);

fig2 = figure('Position',[1351 549 570 413]);
show_target_and_tracker(fig2, target, xy_dpekf_hist,eig_Pdpekf_est_hist,eig_Pdpekf_pred_hist, 'DPEKF');

fig3 = figure('Position',[74 1 570 413]);
show_target_and_tracker_dist(fig3, target, xy_sspekf_hist,eig_Psspekf_est_hist,eig_Psspekf_pred_hist, 'SSPEKF');

fig4 = figure('Position',[1943 38 560 420]);
kfs = {xy_dpekf_hist,xy_sspekf_hist(:,1,:),xy_sspekf_hist(:,2,:),xy_sspekf_hist(:,3,:),xy_sspekf_hist(:,4,:)};
kfs_labels = {'DPEKF','SSPEKF-bs1','SSPEKF-bs2','SSPEKF-bs3','SSPEKF-bs4'};
show_rmse_comparison(fig4,target,kfs, kfs_labels);