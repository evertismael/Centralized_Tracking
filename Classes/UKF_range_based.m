classdef UKF_range_based
   properties
      x_est
      P_est
      x_pred
      P_pred
      eig_P_est
      eig_P_pred
      
      f
      A
      Q
      var_v
      
      h
      R
      bx
      
      alpha
      beta
      kappa
      sig_wc
      sig_wm
   end
   methods
      function obj = UKF_range_based()
          scene = Params.get_scene();
          obj.bx = scene.bx;
          
          % define f as constant velocity model:
          [kf_params, ekf_params, ukf_params] = Params.get_kf_params();
          
          obj.var_v = kf_params.var_v;
          obj.A = @(T) [1 T 0 0;
                        0 1 0 0;
                        0 0 1 T;
                        0 0 0 1];
          obj.f = @(T,x)obj.A(T)*x;
          Qx = @(T,var_v) [0.5*(T^2); T]*var_v*([0.5*(T^2); T].');
          obj.Q = @(T,var_v)[Qx(T,var_v) zeros(2,2); zeros(2,2) Qx(T,var_v)];
          
          % define h for range measurements:
          obj.R = @(deltas_var) diag(deltas_var);
          A1 = [1 0 0 0; 0 0 1 0];
          obj.h = @(x_vect) sqrt(sum((A1*x_vect - obj.bx).^2,1)); 
          
          % get sigma weights
          obj.alpha = ukf_params.alpha;
          obj.beta = ukf_params.beta;
          obj.kappa = ukf_params.kappa;
          n = 4;
          [obj.sig_wm, obj.sig_wc] = sigma_weights(n,obj.alpha, obj.beta, obj.kappa);
          '';
      end
      
      function obj = set_x0(obj,x_0)
          obj.x_est = ([1 0 0 0; 0 0 1 0].')*x_0;
          obj.x_est = obj.x_est + [0,rand(1),0,rand(1)].';
          obj.P_est = 10*eye(4);
          obj.eig_P_est = eig(obj.P_est);
          obj.eig_P_pred = obj.eig_P_est;
          '';
      end
      
      function obj = predict(obj, T)
          obj.x_pred = obj.f(T,obj.x_est);
          obj.P_pred = obj.A(T)*obj.P_est*(obj.A(T).') + obj.Q(T,obj.var_v);
          obj.eig_P_pred = eig(obj.P_pred);
      end
      
      function obj = correct(obj, z_mean, z_var)
          % sigma points
          n = size(obj.x_est,1);
          sig_x = sigma_points(obj.alpha, obj.kappa, obj.x_pred, obj.P_pred);
                    
          % compute z_ut mean:
          sig_z = zeros(size(obj.bx,2),2*n+1);
          for sig_idx = 1:2*n+1
              sig_z(:,sig_idx) = obj.h(sig_x(:,sig_idx));
          end
          z_ut = sum(obj.sig_wm.*sig_z,2);
          
          % compute z_ut var:
          Pz_ut = zeros(size(sig_z,1),size(sig_z,1));
          for sig_idx = 1:2*n+1
             dz = (sig_z(:,sig_idx) - z_ut);
             Pz_ut = Pz_ut + obj.sig_wc(:,sig_idx)*(dz*dz.');
          end
          Pz_ut = Pz_ut + obj.R(z_var);
          '';
          % compute xz_ut_var
          Pxz_ut = zeros(size(sig_x,1),size(sig_z,1));
          for sig_idx = 1:2*n+1
             dx = (sig_x(:,sig_idx) - obj.x_pred);
             dz = (sig_z(:,sig_idx) - z_ut);
             Pxz_ut = Pxz_ut + obj.sig_wc(:,sig_idx)*(dx*dz.');
          end
          
          
          % kalman gain:
          K = Pxz_ut*pinv(Pz_ut);
          
          % correction:
          obj.x_est = obj.x_pred + K*(z_mean - z_ut);
          obj.P_est = obj.P_pred - K*Pz_ut*(K.');
          obj.eig_P_est = eig(obj.P_est);
      end

   end
end