function [mean_xy, var_xy] = sspe_vi_dot_moments(vi_dot)
    scene = Params.get_scene();
    gs = Params.get_grid_search();
    
    lposterior = scene.N_bs*vi_dot;
    posterior_xy = exp(lposterior - max(lposterior,[],[1,2,3]));
    post_norm = trapz(trapz(posterior_xy,1),2)*gs.dx*gs.dy;
    posterior_x = trapz(posterior_xy,2)*gs.dy;
    posterior_y = trapz(posterior_xy,1)*gs.dx;

    % means:
    x_mean_v = trapz(gs.x.*posterior_x,1)*gs.dx./post_norm;
    y_mean_v = trapz(gs.y.*posterior_y,2)*gs.dy./post_norm;

    % variances:
    var_x_v = trapz((gs.x-x_mean_v).^2.*posterior_x,1)*gs.dx./post_norm;
    var_y_v= trapz((gs.y-y_mean_v).^2.*posterior_y,2)*gs.dy./post_norm;
    var_xy_v = trapz((gs.x-x_mean_v).*trapz((gs.y-y_mean_v).*posterior_xy,2),1)*gs.dx*gs.dy./post_norm;
    
    % output:
    mean_xy = zeros(2,scene.N_bs);
    var_xy = zeros(2,2,scene.N_bs);
    for bs_idx=1:scene.N_bs
      mean_xy(:,bs_idx) = [x_mean_v(bs_idx) y_mean_v(bs_idx)].';
      var_xy(:,:,bs_idx) = [var_x_v(bs_idx) var_xy_v(bs_idx); var_xy_v(bs_idx) var_y_v(bs_idx)];
    end
    '';
end
