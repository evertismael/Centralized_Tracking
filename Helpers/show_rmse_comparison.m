function show_rmse_comparison(fig,target,kfs, kfs_labels)
figure(fig);
hold on; grid on;
for idx = 1:size(kfs,2)
    xy_true = target.history([1,3],:);
    xy_hat = kfs{1,idx}([1,3],:);
    rmse = sqrt(sum((xy_true - xy_hat).^2,1));
    plot(target.t_vect,rmse,'DisplayName',kfs_labels{1,idx});
    
end
legend();
end