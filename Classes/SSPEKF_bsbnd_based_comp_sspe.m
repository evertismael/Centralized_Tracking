classdef SSPEKF_bsbnd_based_comp_sspe
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
      function obj = SSPEKF_bsbnd_based_comp_sspe()
          
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
          % SSPE-prior: cnss on prior means and variances:
          vi_means = obj.x_pred_v;
          vi_vars = obj.P_pred_v;
          
          sum_vj_means = zeros(size(vi_means));          
          sum_vj_vars = zeros(size(vi_vars));          
          
          for iter_idx = 1:obj.N_iter
              vj_means=vi_means;
              vj_vars=vi_vars;
              
              for bs_idx=1:scene.N_bs
                  sum_vj_means(:,bs_idx) = sum(vj_means(:,network.A(bs_idx,:)),2);
                  sum_vj_vars(:,:,bs_idx) = sum(vj_vars(:,:,network.A(bs_idx,:)),3);
              end
              
              v_dot_means = obj.x_pred_v + network.beta*(sum_vj_means - network.D.*vi_means);
              v_dot_vars = obj.P_pred_v + network.beta*(sum_vj_vars - reshape(network.D,1,1,scene.N_bs).*vi_vars);
              
              vi_means = vi_means + v_dot_means;
              vi_vars = vi_vars + v_dot_vars;
          end
          
          prior_means = v_dot_means;
          prior_vars = v_dot_vars;
            
          % generate prior:
          lprior_grid = zeros(size(llk_grid));
          for bs_idx = 1:scene.N_bs
              P_pred_inv = pinv(prior_vars(:,:,bs_idx));
              diff_x = gs.x - prior_means(1,bs_idx);
              diff_y = gs.y - prior_means(2,bs_idx); 
              log_prior_bs = P_pred_inv(1,1)*diff_x.^2 +...
                        (P_pred_inv(1,2) +  P_pred_inv(2,1))*diff_x.*diff_y + ...
                         P_pred_inv(2,2)*diff_y.^2;
              
              lprior_grid(:,:,1,bs_idx) = -0.5*log_prior_bs;
          end
                    
          '';
          % compute mean and variance:
          [mean_xy, var_xy] = sspe_vi_dot_moments(v_dot + (1/scene.N_bs)*lprior_grid);
          
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