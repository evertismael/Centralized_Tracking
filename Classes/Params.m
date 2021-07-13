classdef Params
   methods (Static)
       function scene = get_scene()
           persistent bx;
           persistent network_type;
           persistent sim_name;
           
           network_type = 'fully'; % fully / ring / default
           sim_name = 'person_oval'; % car_turn / person_oval
           
           scene = {};
           scene.sim_name = sim_name;
           if strcmp(sim_name,'person_oval')
               bx = [0,0; 50,0; 50,50; 0,50].';
               scene.bx = bx;
               scene.N_bs = size(bx,2);
               % Adj matrix:
               if strcmp(network_type,'fully')
                   % Fully connected netxork:
                   scene.bs_from = [1 1 1, 2 2 2, 3 3 3, 4 4 4];
                   scene.bs_to  =  [2 3 4, 1 3 4, 1 2 4, 1 2 3];
                   scene.A = zeros(scene.N_bs,scene.N_bs);
                   idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
                   scene.A(idx_from_to) = 1;
               elseif strcmp(network_type,'ring')
                   % Ring-Gossip network:
                   scene.bs_from = [1 1, 2 2, 3 3, 4 4];
                   scene.bs_to  =  [2 4, 1 3, 2 4, 1 3];
                   scene.A = zeros(scene.N_bs,scene.N_bs);
                   idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
                   scene.A(idx_from_to) = 1;
               else
                   error('network_type not recognized');
               end
           elseif strcmp(sim_name,'car_turn')
               bx = [30,90; 60,90; 0,60; 30,60; 60,60; 90,60; 0,30; 30,30; 60,30; 90,30; 30,0; 60,0].';
               scene.bx = bx;
               scene.N_bs = size(bx,2);
               % Adj matrix:
               if strcmp(network_type,'fully')
                   % Fully connected netxork:
                   scene.bs_from = [];
                   scene.bs_to  =  [];
                   all_nodes = 1:scene.N_bs;
                   for bs_idx = 1:scene.N_bs
                       scene.bs_from = [scene.bs_from, bs_idx*ones(1,scene.N_bs-1)];
                       scene.bs_to = [scene.bs_to, all_nodes(1:end~=bs_idx)];
                   end
                   scene.A = zeros(scene.N_bs,scene.N_bs);
                   idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
                   scene.A(idx_from_to) = 1;
               elseif strcmp(network_type,'default')
                   % Ring-Gossip network:
                   scene.bs_from = [1 1, 2 2, 3 3, 4 4 4 4, 5 5 5 5, 6 6 , 7 7, 8 8 8 8,  9 9 9 9,   10 10, 11 11, 12 12];
                   scene.bs_to  =  [2 4, 1 5, 7 4, 1 3 5 8, 2 4 6 9, 5 10, 3 8, 4 7 9 11, 5 8 10 12, 6 9, 8 12, 9 11];
                   scene.A = zeros(scene.N_bs,scene.N_bs);
                   idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
                   scene.A(idx_from_to) = 1;
               else
                   error('network_type not recognized');
               end
           else
               error('sim_name not recognized');
           end
       end
       
       function gs = get_grid_search()
           scene = Params.get_scene();
           gs = {};
           gs.x_min = 0; gs.x_max = max(scene.bx(1,:),[],'all'); gs.Nx = 50;
           gs.y_min = 0; gs.y_max = max(scene.bx(2,:),[],'all'); gs.Ny = 50;
           gs.r_min = 0; gs.r_max = 150; gs.Nr = 200;
           
           gs.dx = (gs.x_max - gs.x_min)/gs.Nx;
           gs.dy = (gs.y_max - gs.y_min)/gs.Ny;
           gs.dr = (gs.r_max - gs.r_min)/gs.Nr;
           
           gs.x = (gs.x_min:gs.dx:gs.x_max-gs.dx).';
           gs.y = (gs.y_min:gs.dy:gs.y_max-gs.dy);
           gs.r = (gs.r_min:gs.dr:gs.r_max-gs.dr).';
           
           % 2D R
           gs.delta = sqrt((gs.x - reshape(scene.bx(1,:),1,1,4)).^2 + (gs.y - reshape(scene.bx(2,:),1,1,4)).^2);
           gs.delta = reshape(gs.delta,gs.Nx,gs.Ny,1,4);
       end
       
       function comm = get_communication()
           comm = {};
           comm.c = 2.998e8;
           comm.fc = 2e9;
           comm.B = 40e6;
           comm.N_subcrr = 1024;
           comm.pilot_spacing = 16;
           comm.N_pilot = comm.N_subcrr/comm.pilot_spacing;
           comm.deltaPhi = 2*pi*comm.B/(comm.N_pilot*comm.c);
           comm.Nbps = 2;
           
           comm.pilot_idx_rng = (1:comm.N_pilot);
           comm.phi_rng = comm.deltaPhi*comm.pilot_idx_rng;
           comm.phi_rng = reshape(comm.phi_rng,[1,1,1,1,comm.N_pilot]);
           
           comm.SNR_db = -0;
           comm.noise_type = 'SNR_20m'; % noise_type: SNR_20m/ SNR_center / same
       end
       
       function plt = get_plot()
           scene = Params.get_scene();
           plt = {};
           if strcmp(scene.sim_name,'person_oval')
               plt.xy_axis = [-5 55, -5 55];
           end
       end
       
       function [kf_params, ekf_params, ukf_params] = get_kf_params()
           kf_params = {};
           ekf_params = {};
           ukf_params = {};
           
           sigma_a = (20^2/20)*0.5;
           sigma_v = sqrt(10);
           kf_params.var_a = sigma_a^2;
           kf_params.var_v = sigma_v^2;
           kf_params.N_iter = 1;

           ukf_params.alpha = 0.2; 
           ukf_params.beta = 3-4; %3-n;
           ukf_params.kappa = 2;
           
       end
   end
end