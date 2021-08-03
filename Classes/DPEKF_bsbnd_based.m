classdef DPEKF_bsbnd_based
   properties
      
      x_est
      P_est
      x_pred
      P_pred
      eig_P_est
      eig_P_pred
      
      f
      A
      u
      Q
      var_v
      var_a
      
      h
      H
      bx
      A1

   end
   methods
      function obj = DPEKF_bsbnd_based()
          
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
          
          % define h for range measurements:
          obj.A1 = [1 0; 0 1]; % selector matrix        
     end
      function obj = set_x0(obj,x_0)
          obj.x_est = x_0;
          obj.P_est = 10*eye(2);
          obj.eig_P_est = eig(obj.P_est);
          obj.eig_P_pred = obj.eig_P_est;
          obj.u = 1* rand(2,1);
          '';
     end
      function obj = predict(obj, T)
          obj.x_pred = obj.f(T,obj.x_est) + T*eye(2)*obj.u;
          obj.P_pred = obj.P_est + (T^2)*eye(2)*obj.var_v  + obj.Q(T,obj.var_a);
          obj.eig_P_pred = eig(obj.P_pred);         
          '';
     end
      function obj = correct(obj,bss,T)
          pilot_tx = bss.pilot_tx;
          pilot_rx = bss.pilot_rx;
                   
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          pilot_grid = pilot_tx.*exp(1j*comm.phi_rng.*(gs.delta));
          
          sig_diff = pilot_rx - pilot_grid;
          llk_delta_bs = (-0.5./bss.vars).*(sum(conj(sig_diff).*sig_diff,5));
          
          mem_sum_llk = sum(llk_delta_bs,4);
          
          
          % ADD PREDICTION AS PRIOR:
          P_proy_inv = pinv(obj.P_pred);
          
          diff_x = gs.x - obj.x_pred(1);
          diff_y = gs.y - obj.x_pred(2); 
          
          llk_pred = P_proy_inv(1,1)*diff_x.^2 +...
                    (P_proy_inv(1,2) +  P_proy_inv(2,1))*diff_x.*diff_y + ...
                     P_proy_inv(2,2)*diff_y.^2;
          
          llk_pred = -0.5*llk_pred;
          
          sum_llk = mem_sum_llk + llk_pred;
          
          lk = exp(sum_llk - max(sum_llk,[],'all'));
          lk_norm = trapz(trapz(lk,1),2)*gs.dx*gs.dy;
          lkx = trapz(lk,2)*gs.dy;
          lky = trapz(lk,1)*gs.dx;
          
          % means:
          x_mean = trapz(gs.x.*lkx,1)*gs.dx/lk_norm;
          y_mean = trapz(gs.y.*lky,2)*gs.dy/lk_norm;
          
          
          % variances:
          var_x = trapz((gs.x-x_mean).^2.*lkx,1)*gs.dx/lk_norm;
          var_y = trapz((gs.y-y_mean).^2.*lky,2)*gs.dy/lk_norm;
          var_xy = trapz((gs.x-x_mean).*trapz((gs.y-y_mean).*lk,2),1)*gs.dx*gs.dy/lk_norm;
          
          xy_mean = [x_mean y_mean].';
          xy_var = [var_x var_xy; var_xy var_y];
          
          
          % COMPUTE VELOCITY:
          xy_vel_mean = (xy_mean - obj.x_est)/T;
          
          % correction:
          obj.x_est = xy_mean;
          obj.u = xy_vel_mean;
          obj.P_est = xy_var;
          obj.eig_P_est = eig(obj.P_est);
          '';
      end
   end
end