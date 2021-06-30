function show_target_and_tracker(fig, target, track, eig_P_est_hist, eig_P_pred_hist, tracker_label)
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    plot(target.t_vect,track(1,:),'DisplayName',tracker_label)
    title('x'); ylim([0,Inf]); legend();

    subplot(2,3,4); hold on;
    plot(target.t_vect,target.history(2,:),'DisplayName','true')
    plot(target.t_vect,track(2,:),'DisplayName',tracker_label)
    title('vx');ylim([-25,Inf]); %legend();

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:),'DisplayName','true')
    plot(target.t_vect,track(3,:),'DisplayName',tracker_label)
    title('y'); ylim([0,Inf]); %legend();

    subplot(2,3,5);hold on;
    plot(target.t_vect,target.history(4,:),'DisplayName','true')
    plot(target.t_vect,track(4,:),'DisplayName',tracker_label);
    title('vy');ylim([-25,Inf]); % legend();

    subplot(2,3,6); hold on; grid on;
    plt_params = Params.get_plot();
    axis(plt_params.xy_axis);
    plot(target.history(1,:),target.history(3,:),'--','DisplayName','true','LineWidth',2);
    plot(track(1,:),track(3,:),'DisplayName',tracker_label);
    title('xy-plane'); 
    legend();
    
    subplot(2,3,3); hold on;
    plot(target.t_vect,eig_P_est_hist,'--','DisplayName','est')
    plot(target.t_vect,eig_P_pred_hist,'DisplayName','pred')
    title('eig(P)'); legend();
    
end