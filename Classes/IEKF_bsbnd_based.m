classdef IEKF_bsbnd_based
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
      var_a
      
      h
      H
      R
      bx
      
      A1
      
      % Iterative:
      f_delta_i
      f_grad_d_i
      f_h_i
      N_iter
   end
   methods
      function obj = IEKF_bsbnd_based()
          scene = Params.get_scene();
          obj.bx = scene.bx;
          
          % define f as constant velocity model:
          [kf_params, ~, ~] = Params.get_kf_params();
          obj.var_a = kf_params.var_a;
          obj.var_v = kf_params.var_v;
          obj.A = @(T) [1 T 0 0;
                        0 1 0 0;
                        0 0 1 T;
                        0 0 0 1];
          obj.f = @(T,x)obj.A(T)*x;
          Qx = @(T,var_a) [0.5*(T^2); T]*var_a*([0.5*(T^2); T].');
          obj.Q = @(T,var_a)[Qx(T,var_a) zeros(2,2); zeros(2,2) Qx(T,var_a)];
          
          % define h for range measurements:
          obj.R = @(deltas_var) diag(deltas_var);
          obj.A1 = [1 0 0 0; 0 0 1 0]; % selector matrix
          obj.h = @(x_vect) sqrt(sum((obj.A1*x_vect - obj.bx).^2,1)); 
          
          % linearized version of 'h':
          obj.H = @(x_vect, h_hat) (1./h_hat).*((obj.A1.'*(obj.A1*x_vect - obj.bx)).');
          '';
          
          % Define functions for ITERATIVE EKF:
          obj.f_delta_i = @(x_i) sqrt(sum((x_i - obj.bx).^2,1));
          obj.f_grad_d_i = @(x_i,delta_i) ((x_i - obj.bx)./delta_i).';
          obj.f_h_i = @(x_i,delta_i, phi_rng, pilot_tx) pilot_tx.*exp(1j*phi_rng.*(delta_i));
          obj.N_iter = kf_params.N_iter;
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
          obj.P_pred = obj.A(T)*obj.P_est*(obj.A(T).') + obj.Q(T,obj.var_a);
          obj.eig_P_pred = eig(obj.P_pred);
          
          % check if my derivation is right:
          tmp = obj.A(T);
          Fx = tmp(1:2,1:2);
          Fy = Fx;
          Px = obj.P_est(1:2,1:2);
          Py = obj.P_est(3:4,3:4);
          Pxy = obj.P_est(1:2,3:4);
          Pyx = obj.P_est(3:4,1:2);
          
          Qx = @(T,var_a) [0.5*(T^2); T]*var_a*([0.5*(T^2); T].');
          G = Qx(T,obj.var_a);
          t11 = Fx*Px*Fx.' + G;
          t12 = Fx*Pxy*Fy.';
          t21 = Fy*Pyx*Fx.';
          t22 = Fy*Py*Fy.' + G;
          
          d11 = det(t11);
          d12 = det(t12);
          d21 = det(t21);
          d22 = det(t22);
          obj.P_pred ;
          [d11, d12 d21 d22];
          '';
      end

      function obj = correct(obj, z_mean, z_var, pilot_tx)
          % init values and const inverse matrices:
          x_i = obj.x_pred;
          comm = Params.get_communication();
          phi_rng = squeeze(comm.phi_rng);
          x_i_hist = zeros(4,obj.N_iter);
          
          tmp_c = repmat(z_var.',64,1);
          inv_C = diag(1./tmp_c(:));
          inv_R_1 = 2*[inv_C, zeros(size(inv_C)); zeros(size(inv_C)), inv_C];
          
          for iter_idx = 1:obj.N_iter
              d_i = obj.f_delta_i(obj.A1*x_i);
              grad_d_i = obj.f_grad_d_i(obj.A1*x_i,d_i);
              h_i = obj.f_h_i(obj.A1*x_i,d_i,phi_rng,pilot_tx);
              
              % compute H at x_i
              tmp_h = 1j*phi_rng.*h_i;
              tmp_1 = tmp_h.*(grad_d_i(:,1).');
              tmp_2 = tmp_h.*(grad_d_i(:,2).');
              grad_h = [tmp_1(:), zeros(size(tmp_1(:))),tmp_2(:), zeros(size(tmp_1(:)))];
              H_1 = [real(grad_h); imag(grad_h)];
              
              % y at x_i
              tmp_y = z_mean.' - h_i;
              tmp_y = tmp_y(:);
              y = [real(tmp_y); imag(tmp_y)];
              
              % kalman gain & update
              K = pinv((H_1.')*inv_R_1*H_1 + pinv(obj.P_pred))*(H_1.')*inv_R_1;
              x_i_new = obj.x_pred + K*(y - H_1*(obj.x_pred-x_i));
              
              % save history
              x_i_hist(:,iter_idx) = x_i;
              x_i = x_i_new;
          end  
                   
          P_i = pinv((H_1.')*inv_R_1*H_1 + pinv(obj.P_pred));
            
          % correction:
          obj.x_est = x_i;
          obj.P_est = P_i;
          obj.eig_P_est = eig(obj.P_est);
      end
      
      
      
            
%       function obj = correct_just_hessian(obj, z_mean, z_var, pilot_tx,bss)
%           % init values and const inverse matrices:
%           x_i = obj.x_pred;
%           
%           tmp_c = repmat(z_var.',64,1);
%           inv_C = diag(1./tmp_c(:));
%           inv_R_1 = 2*[inv_C, zeros(size(inv_C)); zeros(size(inv_C)), inv_C];
%           
%           inv_P = pinv(obj.P_pred);
%           
%           comm = Params.get_communication();
%           phi_rng = squeeze(comm.phi_rng);
%           
%           N_iter = 20;
%           x_i_hist = zeros(4,N_iter);
%           for iter_idx = 1:N_iter
%               d_i = obj.f_delta_i(obj.A1*x_i);
%               grad_d_i = obj.f_grad_d_i(obj.A1*x_i,d_i);
%               h_i = obj.f_h_i(obj.A1*x_i,d_i,phi_rng,pilot_tx);
%               % compute H at x_i
%               tmp_h = 1j*phi_rng.*h_i;
%               tmp_1 = tmp_h.*(grad_d_i(:,1).');
%               tmp_2 = tmp_h.*(grad_d_i(:,2).');
%               grad_h = [tmp_1(:), zeros(size(tmp_1(:))),tmp_2(:), zeros(size(tmp_1(:)))];
%               H_1 = [real(grad_h); imag(grad_h)];
%               
%               tmp_y = z_mean.' - h_i;
%               tmp_y = tmp_y(:);
%               y = [real(tmp_y); imag(tmp_y)];
%               
%               % compute gradient:
%               grad = -(H_1.')*inv_R_1*y + inv_P*(x_i - obj.x_pred);
%               hess_inv = pinv((H_1.')*inv_R_1*H_1 + inv_P);
%               
%               x_i_new = x_i - hess_inv*grad;
%               
%               x_i_hist(:,iter_idx) = x_i;
%               x_i = x_i_new;
%           end  
% %           figure;
% %           plot(x_i_hist.');
% %           legend();
% %           
%           P_i = pinv((H_1.')*inv_R_1*H_1 + pinv(obj.P_pred));
%             
%           % correction:
%           obj.x_est = x_i;
%           obj.P_est = P_i;
%           obj.eig_P_est = eig(obj.P_est);
%           
%       end
      
      
   end
end