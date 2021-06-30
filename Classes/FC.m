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
      function [xy_tdoa, varxy_tdoa] = multilateration_tdoa(obj, deltas_mean, deltas_var)
          scene = Params.get_scene();
          gs = Params.get_grid_search();
          tmp_x = gs.x-reshape(obj.bx(1,:),1,1,scene.N_bs);
          tmp_y = gs.y-reshape(obj.bx(2,:),1,1,scene.N_bs);
          
          gsr = sqrt(tmp_x.^2 + tmp_y.^2);
          
          % assume a time synchronized system => t0 = 0;
          % compute position based on gird-search
          deltas_mean = reshape(deltas_mean,1,1,scene.N_bs);
          deltas_var = reshape(deltas_var,1,1,scene.N_bs);
          
          [~,dmin_idx] = min(deltas_mean);
          %dmin_idx
          
          diff_deltas_mean = deltas_mean(1:end~=dmin_idx) - deltas_mean(dmin_idx);
          diff_deltas_var = deltas_var(1:end~=dmin_idx) + deltas_var(dmin_idx);
          
          % diff_grid:
          diff_gsr = gsr(:,:,1:end~=dmin_idx) - gsr(:,:,dmin_idx);
          
          llk = (-0.5./diff_deltas_var).*(diff_gsr-diff_deltas_mean).^2;
          sum_llk = sum(llk,3);
          
          % get pdf:
          lk = exp(sum_llk-max(sum_llk,[],'all'));
          lk_y = trapz(lk,1)*gs.dx;
          lk_x = trapz(lk,2)*gs.dy;
          pdf_norm = trapz(lk_y,2)*gs.dy;
          
          % compute mean and var
          mean_x = (1/pdf_norm)*trapz(gs.x.*lk_x,1)*gs.dx;
          mean_y = (1/pdf_norm)*trapz(gs.y.*lk_y,2)*gs.dy;
          xy_tdoa = [mean_x,mean_y].';
          
          tmp_y = trapz(((gs.x - mean_x)).*lk,1)*gs.dx; % for easy computation
          var_x = (1/pdf_norm)*trapz((gs.x - mean_x).^2.*lk_x,1)*gs.dx;
          var_y = (1/pdf_norm)*trapz((gs.y - mean_y).^2.*lk_y,2)*gs.dy;
          var_xy = (1/pdf_norm)*trapz((gs.y - mean_y).*tmp_y,2)*gs.dy;
          varxy_tdoa = [var_x var_xy; var_xy var_y];  
          
%           figure;
%           imagesc(lk);
%           title('tdoa');
%           ''
      end
      
      % PRIOR: IPE:
      function [prior_mean, prior_var] = prior_toa(obj, xy_mean, xy_var, act_bss)
          % prior mean:
          prior_mean = sqrt(sum((xy_mean - obj.bx(:,act_bss)).^2,1));
          
          % prior variance:
          xy_diff = obj.bx(:,act_bss) - xy_mean;
          tmp = [sum(xy_diff.*xy_var(:,1),1); sum(xy_diff.*xy_var(:,2),1)];
          prior_var = (1./(prior_mean.^2)).*sum(tmp.*xy_diff,1);
          
          % check small variance:
          prior_var(prior_var<1e-10) = 1e-10;
          '';
      end      
   end
end