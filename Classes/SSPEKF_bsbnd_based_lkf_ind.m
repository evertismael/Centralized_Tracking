classdef SSPEKF_bsbnd_based_lkf_ind
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
      function obj = SSPEKF_bsbnd_based_lkf_ind()
          
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
          
          % ---------------------------------------------------------------
          % SSPE: only the measurements (bsbnd signals):
          vi = llk_grid;
          sum_vj = zeros(size(vi));          
          msg = sspe_vi_compress(vi);
          for iter_idx = 1:obj.N_iter
              vj = sspe_vi_decompress(msg);
              for bs_idx=1:scene.N_bs
                  sum_vj(:,:,1,bs_idx) = sum(vj(:,:,:,network.A(bs_idx,:)),4);
              end
              v_dot = llk_grid + network.beta*(sum_vj - reshape(network.D,1,1,1,scene.N_bs).*vi);
              vi = vi + v_dot;
              msg = sspe_vi_compress(vi);
          end
          % ---------------------------------------------------------------
          % compute mean and variance:
          [mean_z, var_z] = sspe_vi_dot_moments(v_dot);
          
          % Direct independent lkf:
          mean_xy = zeros(2,scene.N_bs);
          var_xy = zeros(2,2,scene.N_bs);
          
          for bs_idx = 1:scene.N_bs
            y = mean_z(:,bs_idx) - obj.x_pred_v(:,bs_idx);
            H_i = eye(2);
            S = H_i*obj.P_pred_v(:,:,bs_idx)*(H_i.') + var_z(:,:,bs_idx);
            K = obj.P_pred_v(:,:,bs_idx)*(H_i.')*pinv(S);
            
            mean_xy(:,bs_idx) = obj.x_pred_v(:,bs_idx) + K*y;
            var_xy(:,:,bs_idx) = (eye(2)-K*H_i)*obj.P_pred_v(:,:,bs_idx);
          end
          
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