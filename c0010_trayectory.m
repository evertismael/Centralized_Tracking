clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
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

scene = Params.get_scene();
N_t = size(target.t_vect,2);

deltas_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
for t_idx = 1:N_t
    target.history(:,t_idx);
    
    % noise_type: SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx));
    [~, deltas_mean, deltas_var,act_bss] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
    
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    xy_toa_hist(:,t_idx) = xy_toa;
end

fig1 = figure('Position',[1925 847 560 420]);
show_target_and_xy_toa_in_state_space(fig1,target,xy_toa_hist);

fig2 = figure('Position',[1943 38 560 420]);
show_target_meas_xy_toa_in_meas_space(fig2,target,xy_toa_hist, deltas_hist);
