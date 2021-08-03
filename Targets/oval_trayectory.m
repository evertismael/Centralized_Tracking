function target = oval_trayectory()
x0 = [25,20,5,0].';
target = Mobile(x0);
syms t

% constant velocity model: but we could have a non-linear system here
R = 20;
w0 = sqrt(sum(x0([2,4]).^2))/R;
t1 = 0; t2 = 2*2*pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init));
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

end