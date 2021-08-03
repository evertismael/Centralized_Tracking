function msg = sspe_vi_compress(vi)
    scene = Params.get_scene();
    gs = Params.get_grid_search();
    
    % normalization:
    pdf_vi = exp(vi - max(vi,[],[1,2,3]));
    pdf_norm = trapz(trapz(pdf_vi,1),2)*gs.dx*gs.dy;
    
    % means:
    rng_mean = trapz(trapz(gs.delta.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    tht_mean = trapz(trapz(gs.theta.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    % vars:
    rng_var = trapz(trapz((gs.delta-rng_mean).^2.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    diff_tht = (gs.theta-tht_mean);
    diff_tht(diff_tht >= pi) = diff_tht(diff_tht >= pi) - 2*pi;
    diff_tht(diff_tht <= -pi) = diff_tht(diff_tht <= -pi) + 2*pi;
    tht_var = trapz(trapz(diff_tht.^2.*pdf_vi,1),2)*gs.dx*gs.dy./pdf_norm;
    
    rng_var(rng_var<1e-30) = 1e-30;
    tht_var(tht_var<1e-30) = 1e-30;
    

    msg = {};
    msg.rng_mean = rng_mean;
    msg.tht_mean = tht_mean;
    msg.rng_var = rng_var;
    msg.tht_var = tht_var;
    '';
end