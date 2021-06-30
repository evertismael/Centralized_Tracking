function show_target_meas_xy_toa_in_meas_space(fig, target, xy_toa_hist, deltas_hist)
scene = Params.get_scene();

% proyect target to meas space:
xy_target = target.history([1,3],:);
xy_target = reshape(xy_target,2,1,size(xy_target,2));
target_proj = squeeze(sqrt(sum((xy_target - scene.bx).^2, 1)));

xy_toa = reshape(xy_toa_hist,2,1,size(xy_toa_hist,2));
toa_proj = squeeze(sqrt(sum((xy_toa - scene.bx).^2, 1)));





figure(fig)
    for bs_idx = 1:scene.N_bs
        subplot(2,3,bs_idx); hold on;
        plot(target.t_vect,target_proj(bs_idx,:),'DisplayName','true')
        plot(target.t_vect,toa_proj(bs_idx,:),'DisplayName','xy-toa')
        plot(target.t_vect,deltas_hist(bs_idx,:),'DisplayName','meas-bsbnd')
        title(['BS-',num2str(bs_idx)]); 
        if bs_idx == 1
            legend();
        end
    end
end