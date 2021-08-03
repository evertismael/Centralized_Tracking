classdef Params
   methods (Static)
       function scene = get_scene()
           persistent bx;
           persistent sim_name;
           
           sim_name = 'person_oval'; % car_turn / person_oval
           
           scene = {};
           scene.sim_name = sim_name;
           if strcmp(sim_name,'person_oval')
               bx = [0,0; 50,0; 50,50; 0,50].';
           elseif strcmp(sim_name,'car_turn')
               bx = [30,90; 60,90; 0,60; 30,60; 60,60; 90,60; 0,30; 30,30; 60,30; 90,30; 30,0; 60,0].';
           else
               error('sim_name not recognized');
           end
           
           scene.bx = bx;
           scene.N_bs = size(bx,2);
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
           
           % 2D R:
           gs.delta = sqrt((gs.x - reshape(scene.bx(1,:),1,1,scene.N_bs)).^2 + (gs.y - reshape(scene.bx(2,:),1,1,scene.N_bs)).^2);
           gs.delta = reshape(gs.delta,gs.Nx,gs.Ny,1,scene.N_bs);
           
           % 2D theta:
           gs.theta = atan2((gs.y - reshape(scene.bx(2,:),1,1,scene.N_bs)), (gs.x - reshape(scene.bx(1,:),1,1,scene.N_bs)));
           gs.theta = reshape(gs.theta,gs.Nx,gs.Ny,1,scene.N_bs);
           '';
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
           
           comm.SNR_db = -10;
           comm.noise_type = 'SNR_20m'; % noise_type: SNR_20m/ SNR_center / same
       end
       
       function network = get_network()
           scene = Params.get_scene();
           network = {};
           
           persistent network_type;
           network_type = 'ring'; % fully / ring / default
           
           % Adj matrix:
           if strcmp(network_type,'fully')
               network.A = ones(scene.N_bs) - eye(scene.N_bs);
               '';
           elseif strcmp(network_type,'ring')
               if strcmp(scene.sim_name,'person_oval')
                   network.bs_from = [1 1, 2 2, 3 3, 4 4];
                   network.bs_to  =  [2 4, 1 3, 2 4, 1 3];
               elseif strcmp(scene.sim_name,'car_turn')
                   network.bs_from = [1 1, 2 2, 3 3, 4 4 4 4, 5 5 5 5, 6 6 , 7 7, 8 8 8 8,  9 9 9 9,   10 10, 11 11, 12 12];
                   network.bs_to  =  [2 4, 1 5, 7 4, 1 3 5 8, 2 4 6 9, 5 10, 3 8, 4 7 9 11, 5 8 10 12, 6 9, 8 12, 9 11];
               end
               network.A = zeros(scene.N_bs,scene.N_bs);
               idx_from_to = sub2ind(size(network.A), network.bs_from,network.bs_to);
               network.A(idx_from_to) = 1;
               '';
           else
               error('network_type not recognized');
           end
           network.A = logical(network.A);
           network.D = reshape(sum(network.A,2),1,scene.N_bs);
           network.beta = 0.15;
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
           kf_params.N_iter = 3;

           ukf_params.alpha = 0.2; 
           ukf_params.beta = 3-4; %3-n;
           ukf_params.kappa = 2;
           
       end
   end
end