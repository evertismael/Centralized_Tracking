function show_target_and_tracker_dist(fig, target, tracks, eig_Ps_est_hist, eig_Ps_pred_hist, tracker_label)
N_bs = size(tracks,2);
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    for i=1:N_bs
        plot(target.t_vect,squeeze(tracks(1,i,:)),'DisplayName',[tracker_label,'-',num2str(i)]);
    end
    title('x'); ylim([0,Inf]); legend();

    subplot(2,3,4); hold on;
    plot(target.t_vect,target.history(2,:),'DisplayName','true')
    for i=1:N_bs
        plot(target.t_vect,squeeze(tracks(2,i,:)),'DisplayName',[tracker_label,'-',num2str(i)]);
    end
    title('vx');ylim([-25,Inf]); %legend();

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:),'DisplayName','true')
    for i=1:N_bs
        plot(target.t_vect,squeeze(tracks(3,i,:)),'DisplayName',[tracker_label,'-',num2str(i)]);
    end
    title('y'); ylim([0,Inf]); %legend();

    subplot(2,3,5);hold on;
    plot(target.t_vect,target.history(4,:),'DisplayName','true')
    for i=1:N_bs
        plot(target.t_vect,squeeze(tracks(4,i,:)),'DisplayName',[tracker_label,'-',num2str(i)]);
    end
    title('vy');ylim([-25,Inf]); % legend();

    subplot(2,3,6); hold on; grid on;
    draw_scene();
    plt_params = Params.get_plot();
    axis(plt_params.xy_axis);
    plot(target.history(1,:),target.history(3,:),'--','DisplayName','true','LineWidth',2);
    for i=1:N_bs
        plot(squeeze(tracks(1,i,:)),squeeze(tracks(3,i,:)),'DisplayName',[tracker_label,'-',num2str(i)]);
    end
    title('xy-plane');
    %legend();
    
    subplot(2,3,3); hold on;
    for i=1:N_bs
        plot(target.t_vect,squeeze(eig_Ps_est_hist(:,i,:)),'--','DisplayName',['est-bs',num2str(i)]);
        plot(target.t_vect,squeeze(eig_Ps_pred_hist(:,i,:)),'DisplayName',['pred',num2str(i)]);
    end
    title('eig(P)'); %legend();
    
end