function draw_roads()
    % BSs:
    scene = Params.get_scene();
    scatter(scene.bx(1,:),scene.bx(2,:),'k^');
    text(scene.bx(1,:)+1,scene.bx(2,:)+1,string(1:scene.N_bs));
    
    % Road:
    plot([-5 3 7 7],[43 43 47 55],'k');
    plot([-5 3 7 7],[7 7 3 -5],'k');
    plot([43 43 47 55],[-5 3 7 7],'k');
    plot([43 43 47 55],[55 47 43 43],'k');
    
    plot([25 25],[-5, 0],'k--');
    plot([25 25],[50, 55],'k--');
    
    plot([-5, 0],[25 25],'k--');
    plot([50, 55],[25 25],'k--');
    
    
    % roundabout:
    draw_circle(25,25,17); 

end