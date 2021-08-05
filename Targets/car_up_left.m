function target = car_up_left()
x0 = [30,0,2,15].';
target = Mobile(x0);
syms t

% TRAYECTORY 1: up 
t1 = 0; t2 = 0.9;
a = 12; b = 0.5; c = .3;
law_vx = @(x_init,v_init,t_init,t) v_init(1) +  0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) - b./(1+exp(-a.*(t-c)))*v_init(2);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 2: turning to right to start the roundabout.
w = 1.5;
t1 = t2; t2 = t2 + 1.3;
law_vx = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init) - pi/2)+1.05);
law_vy = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init))+1);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 3: turning to left following the roundabout.
t1 = t2; t2 = t2 + 6.5;
w = 0.51;
phi_x = 0.23*pi;
phi_y = 0.23*pi;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init)+phi_x));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init)+phi_y));
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 4: turning to right to end the roundabout.
t1 = t2; t2 = t2 + 0.8;
w = 3;
law_vx = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init) + pi/20)-2.3);
law_vy = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init) - pi/2)-.5);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

%TRAYECTORY 5: getting away
t1 = t2; t2 = t2 + 1;
a = 16; b = 0.5; c = t1+ .3;
law_vx = @(x_init,v_init,t_init,t) v_init(1) + b./(1+exp(-a.*(t-c)))*v_init(1);
law_vy = @(x_init,v_init,t_init,t) v_init(2)*(1.1^(t-t_init));
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
end