classdef DPEKF_bsbnd_based_avg
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
      A1
      
      x_plot
   end
   methods
      function obj = DPEKF_bsbnd_based_avg()
          
          scene = Params.get_scene();
          obj.bx = scene.bx;
          
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
      function obj = set_x0(obj,x_0)
          scene = Params.get_scene();
          obj.x_est_v = repmat(x_0,1,scene.N_bs);
          obj.P_est_v = repmat(10*eye(2),1,1,scene.N_bs);
          obj.eig_P_est_v = repmat(eig(10*eye(2)),1,scene.N_bs);
          obj.eig_P_pred_v = repmat(eig(10*eye(2)),1,scene.N_bs);
          obj.u_v = 10 * randn(2,scene.N_bs);
          '';
          % extra: x_plot   
          obj.x_plot = [obj.x_est_v(1,:); obj.u_v(1,:); obj.x_est_v(2,:); obj.u_v(2,:)];
          '';
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
      function obj = correct(obj,bss,T)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          scene = Params.get_scene();
          
          pilot_grid = bss.pilot_tx.*exp(1j*comm.phi_rng.*(gs.delta));
          sig_diff = bss.pilot_rx - pilot_grid;
          llk_grid = (-0.5./bss.vars).*(sum(conj(sig_diff).*sig_diff,5));
                                      
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
          
          % add prior to 
          meas_v = llk_grid + (1/scene.N_bs)*lprior_grid;
          
          % DO THIS WITH SSPE
          sspe_vi = mean(meas_v,4);
          '';
         
          % VECTORIZED: mean and variance:
          llk_v = scene.N_bs*sspe_vi;
          lk_v = exp(llk_v - max(llk_v,[],'all'));
          lk_norm_v = trapz(trapz(lk_v,1),2)*gs.dx*gs.dy;
          lkx_v = trapz(lk_v,2)*gs.dy;
          lky_v = trapz(lk_v,1)*gs.dx;
          
          % means:
          x_mean_v = trapz(gs.x.*lkx_v,1)*gs.dx/lk_norm_v;
          y_mean_v = trapz(gs.y.*lky_v,2)*gs.dy/lk_norm_v;
          
          % variances:
          var_x_v = trapz((gs.x-x_mean_v).^2.*lkx_v,1)*gs.dx/lk_norm_v;
          var_y_v= trapz((gs.y-y_mean_v).^2.*lky_v,2)*gs.dy/lk_norm_v;
          var_xy_v = trapz((gs.x-x_mean_v).*trapz((gs.y-y_mean_v).*lk_v,2),1)*gs.dx*gs.dy/lk_norm_v;
          
          xy_mean_v = [x_mean_v y_mean_v].';
          xy_var_v = [var_x_v var_xy_v; var_xy_v var_y_v];
          
          % VECTORIZED: correction:
          xy_vel_mean_v = (xy_mean_v - obj.x_est_v)/T;
          obj.x_est_v = repmat(xy_mean_v,1,scene.N_bs);
          obj.u_v = xy_vel_mean_v;
          obj.P_est_v = repmat(xy_var_v,1,1,scene.N_bs);
          obj.eig_P_est_v = repmat(eig(xy_var_v),1,scene.N_bs);
          
          % extra: x_plot   
          obj.x_plot = [obj.x_est_v(1,:); obj.u_v(1,:); obj.x_est_v(2,:); obj.u_v(2,:)];
          ''; 
      end
   end
end