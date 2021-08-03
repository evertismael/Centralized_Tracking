classdef SSPEKF_bsbnd_based
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
      function obj = SSPEKF_bsbnd_based()
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
          obj.u_v = 10* rand(2,scene.N_bs);
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
          
          % generate prior-grid
          lprior_grid = zeros(size(llk_grid));
          for bs_idx = 1:scene.N_bs
              P_pred_inv = pinv(obj.P_pred_v(:,:,bs_idx));
              diff_x = gs.x - obj.x_pred_v(1,bs_idx);
              diff_y = gs.y - obj.x_pred_v(2,bs_idx); 
              log_prior_bs = P_pred_inv(1,1)*diff_x.^2 +...
                        (P_pred_inv(1,2) +  P_pred_inv(2,1))*diff_x.*diff_y + ...
                         P_pred_inv(2,2)*diff_y.^2;
              
              lprior_grid(:,:,1,bs_idx) = -0.5*log_prior_bs;
          end
          meas_v = llk_grid + (1/scene.N_bs)*lprior_grid;
          
          % SSPE over meas_v (meas + prior):
          vi = meas_v;
          sum_vj = zeros(size(vi));
          for iter_idx = 1:obj.N_iter
              vj = vi;
              for bs_idx=1:scene.N_bs
                  sum_vj(:,:,1,bs_idx) = sum(vj(:,:,:,1:end~=bs_idx),4);
              end
              v_dot = meas_v + network.beta*(sum_vj - scene.N_bs*vi);
              vi = vi + v_dot;
          end
          [xy_mean_v, xy_var_v] = sspe_vi_dot_moments(v_dot);
          '';
          
          % correction:
          xy_vel_mean_v = (xy_mean_v - obj.x_est_v)/T;
          obj.x_est_v = xy_mean_v;
          obj.u_v = xy_vel_mean_v;
          obj.P_est_v = xy_var_v;
          for bs_idx=1:scene.N_bs
            obj.eig_P_est_v(:,bs_idx) = eig(xy_var_v(:,:,bs_idx));
          end
          
          % extra: x_plot
          obj.x_plot = [obj.x_est_v(1,:); obj.u_v(1,:); obj.x_est_v(2,:); obj.u_v(2,:)];
          ''; 
      end
   end
end