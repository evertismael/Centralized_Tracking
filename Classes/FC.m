classdef FC
   properties
      bx
      sigma
   end
   methods
      function obj = FC()
          scene = Params.get_scene();
          obj.bx = scene.bx;
      end
      
      % MULTILATERATION:
      function [xy_toa, varxy_toa] = multilateration_toa(obj, deltas_mean, deltas_var, act_bss)
          gs = Params.get_grid_search();
          tmp_x = gs.x-reshape(obj.bx(1,act_bss),1,1,size(act_bss,2));
          tmp_y = gs.y-reshape(obj.bx(2,act_bss),1,1,size(act_bss,2));
          
          gsr = sqrt(tmp_x.^2 + tmp_y.^2);
          
          % assume a time synchronized system => t0 = 0;
          % compute position based on gird-search
          deltas_mean = reshape(deltas_mean,1,1,size(act_bss,2));
          deltas_var = reshape(deltas_var,1,1,size(act_bss,2));
          
          llk = (-0.5./deltas_var).*(gsr-deltas_mean).^2;
          sum_llk = sum(llk,3);
          
          % get pdf:
          lk = exp(sum_llk-max(sum_llk,[],'all'));
          lk_y = trapz(lk,1)*gs.dx;
          lk_x = trapz(lk,2)*gs.dy;
          pdf_norm = trapz(lk_y,2)*gs.dy;
          
          % compute mean and var
          mean_x = (1/pdf_norm)*trapz(gs.x.*lk_x,1)*gs.dx;
          mean_y = (1/pdf_norm)*trapz(gs.y.*lk_y,2)*gs.dy;
          xy_toa = [mean_x,mean_y].';
          
          tmp_y = trapz(((gs.x - mean_x)).*lk,1)*gs.dx; % for easy computation
          var_x = (1/pdf_norm)*trapz((gs.x - mean_x).^2.*lk_x,1)*gs.dx;
          var_y = (1/pdf_norm)*trapz((gs.y - mean_y).^2.*lk_y,2)*gs.dy;
          var_xy = (1/pdf_norm)*trapz((gs.y - mean_y).*tmp_y,2)*gs.dy;
          varxy_toa = [var_x var_xy; var_xy var_y];  
          
          % check for small variances:
          varxy_toa(varxy_toa < gs.dx^2) = gs.dx^2;
          
%           figure;
%           imagesc(lk);
%           title('toa')
           '';
      end
      
      
      function [xy_mean, xy_var] = dpe(obj,pilot_rx, vars, pilot_tx)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          pilot_grid = pilot_tx.*exp(1j*comm.phi_rng.*(gs.delta));
          
          sig_diff = pilot_rx - pilot_grid;
          llk_delta_bs = (-0.5./vars).*(sum(conj(sig_diff).*sig_diff,5));
          
          sum_llk = sum(llk_delta_bs,4);
          
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
      end
   end
end