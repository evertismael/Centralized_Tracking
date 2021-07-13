function show_target_and_xy_toa_in_state_space(fig, target, xy_toa_hist)
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    plot(target.t_vect,xy_toa_hist(1,:),'DisplayName','multilat')
    %scatter(target.t_vect,xy_toa_hist(1,:),'+','MarkerEdgeAlpha',0.2,'DisplayName','multilat')
    title('x'); ylim([0,Inf]); legend();

    subplot(2,3,4);
    plot(target.t_vect,target.history(2,:))
    title('vx'); ylim([-20,Inf]);

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:));
    plot(target.t_vect,xy_toa_hist(2,:),'DisplayName','multilat')
    %scatter(target.t_vect,xy_toa_hist(2,:),'+','MarkerEdgeAlpha',0.2,'DisplayName','multilat')
    title('y'); ylim([0,Inf]);

    subplot(2,3,5);
    plot(target.t_vect,target.history(4,:))
    title('vy'); ylim([-20,Inf]);

    subplot(2,3,6);
    hold on;
    plot(target.history(1,:),target.history(3,:),'.b');
    scatter(xy_toa_hist(1,:),xy_toa_hist(2,:),'+','MarkerEdgeAlpha',0.7,'DisplayName','multilat')
    plot(target.history(1,1),target.history(3,1),'or'); % initial point
    
    title('xy-plane'); grid on;
    plt_params = Params.get_plot();
    axis(plt_params.xy_axis);
end