classdef SSPEKF_bsbnd_based_opt1
   properties
      x_est_v
      P_est_v
      x_pred_v
      P_pred_v
      eig_P_est_v
      eig_P_pred_v
      
      f
      A
      u_v
      Q
      var_v
      var_a
      
      h
      H
      bx    
      
      N_iter
      x_plot
   end
   methods
      function obj = SSPEKF_bsbnd_based_opt1()
          scene = Params.get_scene();
          kf_params = Params.get_kf_params();
          obj.bx = scene.bx;
          obj.N_iter = kf_params.N_iter;
          
          % define f as constant velocity model:
          [kf_params, ~, ~] = Params.get_kf_params();
          obj.var_v = kf_params.var_v;
          obj.var_a = kf_params.var_a;
          obj.A = @(T) [1 0;
                        0 1];
          obj.f = @(T,x)obj.A(T)*x;
          Qx = @(T,var_a) (0.5*(T^2))*var_a*(0.5*(T^2));
          obj.Q = @(T,var_a)[Qx(T,var_a) 0; 0 Qx(T,var_a)];
     end
      function obj = set_x0(obj,x_0_v)
          scene = Params.get_scene();
          obj.x_est_v = x_0_v;
          obj.P_est_v = repmat(10*eye(2),1,1,scene.N_bs);
          obj.eig_P_est_v = repmat(eig(10*eye(2)),1,scene.N_bs);
          obj.eig_P_pred_v = repmat(eig(10*eye(2)),1,scene.N_bs);
          obj.u_v = 10 * randn(2,scene.N_bs);
          '';
          obj.x_plot = [obj.x_est_v(1,:); obj.u_v(1,:); obj.x_est_v(2,:); obj.u_v(2,:)];
     end
      function obj = predict(obj, T)
          scene = Params.get_scene();
          for bs_idx=1:scene.N_bs
              obj.x_pred_v(:,bs_idx) = obj.f(T,obj.x_est_v(:,bs_idx)) + T*eye(2)*obj.u_v(:,bs_idx);
              obj.P_pred_v(:,:,bs_idx) = obj.P_est_v(:,:,bs_idx) + (T^2)*eye(2)*obj.var_v  + obj.Q(T,obj.var_a);
              obj.eig_P_pred_v(:,bs_idx) = eig(obj.P_pred_v(:,:,bs_idx)); 
          end
          '';
     end
      function obj = correct(obj,bss,T,t_idx)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          scene = Params.get_scene();
          network = Params.get_network();
          
          % generate loglikelihood-grid:
          pilot_grid = bss.pilot_tx.*exp(1j*comm.phi_rng.*(gs.delta));
          sig_diff = bss.pilot_rx - pilot_grid;
          llk_grid = (-0.5./bss.vars).*(sum(conj(sig_diff).*sig_diff,5));
                    
          % generate prior:
          lprior_grid = zeros(size(llk_grid));
          for bs_idx=1:scene.N_bs
              P_proy_inv = pinv(obj.P_pred_v(:,:,bs_idx));
              diff_x = gs.x - obj.x_pred_v(1,bs_idx);
              diff_y = gs.y - obj.x_pred_v(2,bs_idx); 
              llk_pred = P_proy_inv(1,1)*diff_x.^2 +...
                        (P_proy_inv(1,2) +  P_proy_inv(2,1))*diff_x.*diff_y + ...
                         P_proy_inv(2,2)*diff_y.^2;
              lprior_grid(:,:,1,bs_idx) = -0.5*llk_pred;
          end
                    
          meas_v = llk_grid + (1/scene.N_bs)*lprior_grid;
          
          % SSPE: only the measurements (bsbnd signals):
          vi = meas_v;
          sum_vj = zeros(size(vi));          
          
          msg = compress_vi(vi);
          for iter_idx = 1:obj.N_iter
              vj = decompress_vi(msg);
              %vj=vi;
              for bs_idx=1:scene.N_bs
                  sum_vj(:,:,1,bs_idx) = sum(vj(:,:,:,network.A(bs_idx,:)),4);
              end
              v_dot = llk_grid + network.beta*(sum_vj - reshape(network.D,1,1,1,scene.N_bs).*vi);
              vi = vi + v_dot;
              msg = compress_vi(vi);
          end
          
          % compute mean and variance:
          [mean_xy, var_xy] = compute_mean_var(v_dot);
          
          %----------------------------------------------------------------
          % correction:
          xy_vel_mean_v = (mean_xy - obj.x_est_v)/T;
          obj.x_est_v = mean_xy;
          obj.u_v = xy_vel_mean_v;
          obj.P_est_v = var_xy;
          for bs_idx=1:scene.N_bs
            obj.eig_P_est_v(:,bs_idx) = eig(var_xy(:,:,bs_idx));
          end
          
          % extra: x_plot
                  
          obj.x_plot = [obj.x_est_v(1,:); obj.u_v(1,:); obj.x_est_v(2,:); obj.u_v(2,:)];
          ''; 
      end
   end
end

% computes final estimations from delta_vi:
function [mean_xy, var_xy] = compute_mean_var(vi_dot)
    scene = Params.get_scene();
    gs = Params.get_grid_search();
    
    lposterior = scene.N_bs*vi_dot;
    posterior_xy = exp(lposterior - max(lposterior,[],[1,2,3]));
    post_norm = trapz(trapz(posterior_xy,1),2)*gs.dx*gs.dy;
    posterior_x = trapz(posterior_xy,2)*gs.dy;
    posterior_y = trapz(posterior_xy,1)*gs.dx;

    % means:
    x_mean_v = trapz(gs.x.*posterior_x,1)*gs.dx./post_norm;
    y_mean_v = trapz(gs.y.*posterior_y,2)*gs.dy./post_norm;

    % variances:
    var_x_v = trapz((gs.x-x_mean_v).^2.*posterior_x,1)*gs.dx./post_norm;
    var_y_v= trapz((gs.y-y_mean_v).^2.*posterior_y,2)*gs.dy./post_norm;
    var_xy_v = trapz((gs.x-x_mean_v).*trapz((gs.y-y_mean_v).*posterior_xy,2),1)*gs.dx*gs.dy./post_norm;
    
    % output:
    mean_xy = zeros(2,scene.N_bs);
    var_xy = zeros(2,2,scene.N_bs);
    for bs_idx=1:scene.N_bs
      mean_xy(:,bs_idx) = [x_mean_v(bs_idx) y_mean_v(bs_idx)].';
      var_xy(:,:,bs_idx) = [var_x_v(bs_idx) var_xy_v(bs_idx); var_xy_v(bs_idx) var_y_v(bs_idx)];
    end
    '';
end

% compress vi to compute self-synchronization:
function msg = compress_vi(vi)
    scene = Params.get_scene();
    gs = Params.get_grid_search();
    
    % normalization:
    pdf_vi = exp(vi - max(vi,[],[1,2,3]));
    pdf_norm = trapz(trapz(pdf_vi,1),2)*gs.dx*gs.dy;
    
    % means:
    rng_mean = trapz(trapz(gs.delta.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    tht_mean = trapz(trapz(gs.theta.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    % vars:
    rng_var = trapz(trapz((gs.delta-rng_mean).^2.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    diff_tht = (gs.theta-tht_mean);
    diff_tht(diff_tht >= pi) = diff_tht(diff_tht >= pi) - 2*pi;
    diff_tht(diff_tht <= -pi) = diff_tht(diff_tht <= -pi) + 2*pi;
    tht_var = trapz(trapz(diff_tht.^2.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    rng_var(rng_var<1e-30) = 1e-30;
    tht_var(tht_var<1e-30) = 1e-30;
    

    msg = {};
    msg.rng_mean = rng_mean;
    msg.tht_mean = tht_mean;
    msg.rng_var = rng_var;
    msg.tht_var = tht_var;
    '';
end
function vj = decompress_vi(msg)
    gs = Params.get_grid_search();
    
    diff_tht = (gs.theta-msg.tht_mean);
    diff_tht(diff_tht >= pi) = diff_tht(diff_tht >= pi) - 2*pi;
    diff_tht(diff_tht <= -pi) = diff_tht(diff_tht <= -pi) + 2*pi;
    
    vj = -0.5*(1./msg.rng_var).*(gs.delta-msg.rng_mean).^2 + ...
         -0.5*(1./msg.tht_var).*(diff_tht).^2;
end