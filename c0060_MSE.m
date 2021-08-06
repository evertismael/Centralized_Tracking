clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
%rng(2)
% -------------------------------------------------------------------------
% 1.- Define & Generate target trajetories
% -------------------------------------------------------------------------
dt = (1e-1);
%target = oval_trayectory();
target = car_up_left_v2();
target = target.gen_trayectory(dt);

% params and hist
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_sim = 30;

xy_bsbnd_hist = zeros(2,N_t,N_sim);
xy_true_hist = zeros(2,N_t,N_sim);

xy_dpekf_hist = zeros(4,N_t,N_sim);
eig_Pdpekf_est_hist = zeros(2,N_t,N_sim);
eig_Pdpekf_pred_hist = zeros(2,N_t,N_sim);

xy_sspekf_hist = zeros(4,scene.N_bs,N_t,N_sim);
eig_Psspekf_est_hist = zeros(2,scene.N_bs,N_t,N_sim);
eig_Psspekf_pred_hist = zeros(2,scene.N_bs,N_t,N_sim);

xy_sspekf_lkf_hist = zeros(4,scene.N_bs,N_t,N_sim);
eig_Psspekf_lkf_est_hist = zeros(2,scene.N_bs,N_t,N_sim);
eig_Psspekf_lkf_pred_hist = zeros(2,scene.N_bs,N_t,N_sim);

f_rmse = figure();
f_rmsid = figure();
for sim_idx = 1:N_sim
    % -------------------------------------------------------------------------
    % 3.- Define BSs and get ToA measurements from BSs:
    % -------------------------------------------------------------------------
    bss = BSs();
    bss = bss.gen_pilot_tx();
    fc = FC();

    % trackers:
    dpekf = DPEKF_bsbnd_based();
    sspekf = SSPEKF_bsbnd_based_comp_sspe();
    %sspekf = SSPEKF_bsbnd_based_comp_nhop();
    sspekf_lkf = SSPEKF_bsbnd_based_lkf_ind();
    %sspekf_lkf = SSPEKF_bsbnd_based_lkf_all();


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
            sspekf_lkf = sspekf_lkf.set_x0(xy_0_v);
        else
            dpekf = dpekf.predict(dt);
            dpekf = dpekf.correct(bss,dt);

            sspekf = sspekf.predict(dt);
            sspekf = sspekf.correct(bss,dt,t_idx);

            sspekf_lkf = sspekf_lkf.predict(dt);
            sspekf_lkf = sspekf_lkf.correct(bss,dt,t_idx);
        end

        % ______________________________________________
        % save for history:
        % bsbnd:
        xy_true_hist(:,t_idx,sim_idx) = xy_true;
        xy_bsbnd_hist(:,t_idx,sim_idx) = xy_bsbnd;
                
        % dpekf
        xy_dpekf_hist(:,t_idx,sim_idx) = [dpekf.x_est(1) dpekf.u(1) dpekf.x_est(2) dpekf.u(2)].';
        eig_Pdpekf_est_hist(:,t_idx,sim_idx) = dpekf.eig_P_est;
        eig_Pdpekf_pred_hist(:,t_idx,sim_idx) = dpekf.eig_P_pred;

        % sspekf
        xy_sspekf_hist(:,:,t_idx,sim_idx) = sspekf.x_plot;
        eig_Psspekf_est_hist(:,:,t_idx,sim_idx) = sspekf.eig_P_est_v;
        eig_Psspekf_pred_hist(:,:,t_idx,sim_idx) = sspekf.eig_P_pred_v;

        % sspekf_lkf
        xy_sspekf_lkf_hist(:,:,t_idx,sim_idx) = sspekf_lkf.x_plot;
        eig_Psspekf_lkf_est_hist(:,:,t_idx,sim_idx) = sspekf_lkf.eig_P_est_v;
        eig_Psspekf_lkf_pred_hist(:,:,t_idx,sim_idx) = sspekf_lkf.eig_P_pred_v;
    end
    '';
    if mod(sim_idx,3)==0
        % RMSE
        [rmse_bsbnd, rmse_dpekf, rmse_sspekf, rmse_sspekf_lkf] = mse_curves(...
            xy_true_hist(:,:,1:sim_idx),xy_bsbnd_hist(:,:,1:sim_idx),...
            xy_dpekf_hist(:,:,1:sim_idx),xy_sspekf_hist(:,:,:,1:sim_idx),...
            xy_sspekf_lkf_hist(:,:,:,1:sim_idx));
        figure(f_rmse);
        clf();
        hold on; grid;
        plot(target.t_vect, rmse_bsbnd, 'DisplayName','bsbnd');
        plot(target.t_vect, rmse_dpekf, 'DisplayName','dpekf');
        plot(target.t_vect, squeeze(rmse_sspekf), 'DisplayName','sspekf');
        plot(target.t_vect, squeeze(rmse_sspekf_lkf),'--', 'DisplayName','sspekf-lkf');
        legend();
        title(num2str(sim_idx));
        pause(1);
        '';
        
        % RMSID
        rmsid_sspekf = rmsid_curves(xy_sspekf_hist(:,:,:,1:sim_idx));
        rmsid_sspekf_lkf = rmsid_curves(xy_sspekf_lkf_hist(:,:,:,1:sim_idx));
        
        figure(f_rmsid);
        clf();
        hold on; grid;
        plot(target.t_vect, squeeze(rmsid_sspekf), 'DisplayName','sspekf');
        plot(target.t_vect, squeeze(rmsid_sspekf_lkf),'--', 'DisplayName','sspekf-lkf');
        legend();
        title(num2str(sim_idx));
        pause(1);
        '';
    end
end
%%

function [rmse_bsbnd, rmse_dpekf, rmse_sspekf, rmse_sspekf_lkf] = mse_curves(...
    xy_true_hist,xy_bsbnd_hist,xy_dpekf_hist,xy_sspekf_hist,xy_sspekf_lkf_hist)

    % MSE bsbnd:
    tmp = xy_bsbnd_hist - xy_true_hist;
    rmse_bsbnd = sqrt(mean(sum(tmp.^2,1),3));
    
    % MSE dpekf:
    tmp = xy_dpekf_hist([1,3],:,:) - xy_true_hist;
    rmse_dpekf = sqrt(mean(sum(tmp.^2,1),3));
    
    % MSE sspekf:
    Nsim = size(xy_sspekf_hist,4);
    tmp = xy_sspekf_hist([1,3],:,:,:) - reshape(xy_true_hist,2,1,58,Nsim);
    rmse_sspekf = sqrt(mean(sum(tmp.^2,1),4));
    
    % MSE sspekf:
    tmp = xy_sspekf_lkf_hist([1,3],:,:,:) - reshape(xy_true_hist,2,1,58,Nsim);
    rmse_sspekf_lkf = sqrt(mean(sum(tmp.^2,1),4)); 
'';
end

function  rmsid = rmsid_curves(xy_hist)
    '';
    scene = Params.get_scene();
    rmsid = zeros(scene.N_bs,size(xy_hist,3));
    for bs_idx=1:scene.N_bs
        tmp = xy_hist([1,3],1:end~=bs_idx,:,:) - xy_hist([1,3],bs_idx,:,:);
        rmsid(bs_idx,:) = sqrt(mean(mean(sum(tmp.^2,1),4),2));
    end
end