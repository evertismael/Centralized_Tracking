function vj = sspe_vi_decompress(msg)
    gs = Params.get_grid_search();
    
    diff_tht = (gs.theta-msg.tht_mean);
    diff_tht(diff_tht >= pi) = diff_tht(diff_tht >= pi) - 2*pi;
    diff_tht(diff_tht <= -pi) = diff_tht(diff_tht <= -pi) + 2*pi;
    
    vj = -0.5*(1./msg.rng_var).*(gs.delta-msg.rng_mean).^2 + ...
         -0.5*(1./msg.tht_var).*(diff_tht).^2;
end