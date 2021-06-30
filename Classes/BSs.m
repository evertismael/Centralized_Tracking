classdef BSs
   properties
      bx              % positions of BSs
      sigma           % ToA noise of furthest
      SNR_db          % Bsbnd noise of furthest to center position
      c = 2.998e8;    % speed of light
      
      pilot_tx
      pilot_rx
      vars
      active_bs       % BS that are active: detected r 
      pilot_grid_r
      
   end
   methods
      function obj = BSs()
          scene = Params.get_scene();
          comm = Params.get_communication();
          obj.bx = scene.bx;
          obj.SNR_db = comm.SNR_db;
      end
      function obj = gen_pilot_tx(obj)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          % create Tx signal:
          bitstream = randi(2, comm.Nbps * comm.N_pilot, 1) - 1;
          obj.pilot_tx = mapping(bitstream,comm.Nbps,'qam');
          obj.pilot_tx = reshape(obj.pilot_tx,1,1,1,1,comm.N_pilot);
          
          % create signal grid:
          obj.pilot_grid_r = obj.pilot_tx.*exp(1j*comm.phi_rng.*(gs.r));
      end
      function obj = channel_propagation(obj,x_target)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          
          noise_type = comm.noise_type;
          
          x_tx = [x_target(1),x_target(3)].';
          d = sqrt(sum((x_tx - obj.bx).^2,1));
          delta0 = 0; % range offset.
          delta = d + delta0;
          
          % Propagation Channel - Rx signal:
          delta = reshape(delta,1,1,1,length(delta));
          obj.pilot_rx = obj.pilot_tx.*exp(1j*comm.phi_rng.*(delta));
          
          var_n = 10^(-obj.SNR_db/10);
          if strcmp(noise_type,'same')
              w = sqrt(var_n/2)*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));    
              obj.vars = var_n*ones(1,1,1,size(w,4));
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          elseif strcmp(noise_type,'SNR_center')
              % obj.SNR_db: SNR from center to furthest BS
              % distance to center:
              xy_c = [(gs.x_max-gs.x_min)/2, (gs.y_max-gs.y_min)/2].';
              d_c = sqrt(sum((xy_c - obj.bx).^2,1));
              d_c = max(d_c,[],'all');
              % compute SNRs of BSs:
              SNR_lin = 10^(obj.SNR_db/10);
              SNR_i = SNR_lin.*((d_c./d).^2);
              obj.vars = reshape(1./SNR_i,1,1,1,size(d,2));
              w = sqrt(obj.vars/2).*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          elseif strcmp(noise_type,'SNR_20m')
              % obj.SNR_db: SNR when Tx is 20m away from the BS
              d_c = 20;
              % compute SNRs of BSs:
              SNR_lin = 10^(obj.SNR_db/10);
              SNR_i = SNR_lin.*((d_c./d).^2);
              obj.vars = reshape(1./SNR_i,1,1,1,size(d,2));
              w = sqrt(obj.vars/2).*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));
              
              SNR_i_db = 10*log10(SNR_i);
              obj.active_bs = SNR_i_db > -15; % it detects if SNR bigger than 15 db
          else
              w = zeros(size(obj.pilot_rx));
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          end
          obj.pilot_rx = obj.pilot_rx + w;
          
      end
      function [toas, deltas_mean, deltas_var, act_bss] = compute_bsbnd_toa(obj)
          
          gs = Params.get_grid_search();
          
          act_bss = find(obj.active_bs == 1); % active bss
          % compute ToA:
          sig_diff = obj.pilot_rx(:,:,:,act_bss,:) - obj.pilot_grid_r;
          llk_delta_bs = (-0.5./obj.vars(:,:,:,act_bss)).*(sum(conj(sig_diff).*sig_diff,5));
          
          % toa_pdf:
          lk = exp(llk_delta_bs - max(llk_delta_bs,[],[1,2,3]));
          llk_norm = trapz(lk,1)*gs.dr;
          
          deltas_mean = (1/llk_norm).*trapz(gs.r.*lk,1)*gs.dr;
          deltas_var = (1/llk_norm).*trapz(((gs.r-deltas_mean).^2).*lk,1)*gs.dr;
          % check min var:
          deltas_var(deltas_var<gs.dr^2) = gs.dr^2;
          
          % output:
          deltas_mean = squeeze(deltas_mean);
          deltas_var = squeeze(deltas_var);
          toas = deltas_mean./obj.c;
          '';
      end
   end
end