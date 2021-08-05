function target = car_up_left_v2()
x0 = [30,0,2,15].';
target = Mobile(x0);
syms t

% TRAYECTORY 1: turning to right to start the roundabout.
t1 = 0; t2 = 0.32;
R = 4;
a = 0;
h = 35;
k = 3;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(sin((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(cos((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_x = @(x_init,v_init,t_init,t) h - R*(cos((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law_y = @(x_init,v_init,t_init,t) k + R*(sin((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 2: turning to left following the roundabout.
t1 = t2; t2 = t2 + 5.1;
R = 20;
a = 1.5;
h = 25;
k = 25;
law_vx = @(x_init,v_init,t_init,t) -sqrt(sum(v_init.^2,1))*(sin((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_vy = @(x_init,v_init,t_init,t)  sqrt(sum(v_init.^2,1))*(cos((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_x = @(x_init,v_init,t_init,t) h + R*(cos((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law_y = @(x_init,v_init,t_init,t) k + R*(sin((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);


% TRAYECTORY 3: turning to right to start the roundabout.
t1 = t2; t2 = t2 + 0.4;
R = 4;
a = 0.7;
h = 3;
k = 35;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(sin((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(cos((sqrt(sum(v_init.^2,1))/R)*((t-t_init)-a)));
law_x = @(x_init,v_init,t_init,t) h - R*(cos((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law_y = @(x_init,v_init,t_init,t) k + R*(sin((sqrt(sum(v_init.^2,1))/R)*(t-t_init - a)));
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

end