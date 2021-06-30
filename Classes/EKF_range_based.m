classdef EKF_range_based
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
      H
      R
      bx
   end
   methods
      function obj = EKF_range_based()
          scene = Params.get_scene();
          obj.bx = scene.bx;
          
          % define f as constant velocity model:
          kf_params = Params.get_kf_params();
          obj.var_v = kf_params.var_v;
          obj.A = @(T) [1 T 0 0;
                        0 1 0 0;
                        0 0 1 T;
                        0 0 0 1];
          obj.f = @(T,x)obj.A(T)*x;
          Qx = @(T,var_v) [0.5*T^2; T]*var_v*([0.5*T^2; T].');
          obj.Q = @(T,var_v)[Qx(T,var_v) zeros(2,2); zeros(2,2) Qx(T,var_v)];
          
          % define h for range measurements:
          obj.R = @(deltas_var) diag(deltas_var);
          A1 = [1 0 0 0; 0 0 1 0]; % selector matrix
          obj.h = @(x_vect) sqrt(sum((A1*x_vect - obj.bx).^2,1)); 
          
          % linearized version of 'h':
          obj.H = @(x_vect, h_hat) (1./h_hat).*((A1.'*(A1*x_vect - obj.bx)).');
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
          % z_mean, z_var = ranges.
          z_hat = (obj.h(obj.x_pred)).';

          Hk = obj.H(obj.x_pred,z_hat);
                    
          % inovation:
          y = z_mean - z_hat;
          S = Hk*obj.P_pred*(Hk.') + obj.R(z_var);
                    
          % kalman gain:
          K = obj.P_pred*(Hk.')*pinv(S);
          
          % correction:
          obj.x_est = obj.x_pred + K*y;
          obj.P_est = obj.P_pred - K*S*(K.');
          obj.eig_P_est = eig(obj.P_est);
      end
      
   end
end