classdef Mobile
   properties
      x0
      history
      t_vect
      moves_inter
      moves_law
   end
   methods
      function obj = Mobile(x0)
          obj.x0 = x0;
          obj.moves_inter = double.empty(2, 0);
          obj.moves_law = repmat(Law(0,0,0,0),1,0);
      end
      
      function obj = add_trayectory(obj,t1,t2,law)
          obj.moves_inter(:,end+1) = [t1,t2].';
          obj.moves_law(end+1) = law;
      end
      
      function obj = gen_trayectory(obj,dt)
          obj.t_vect = min(obj.moves_inter,[],'all'):dt:max(obj.moves_inter,[],'all')-dt;
          obj.history = zeros(4,size(obj.t_vect,2));
          
          % generate by moves_law
          for law_idx = 1:size(obj.moves_inter,2)
             t0 = obj.moves_inter(1,law_idx);
             t_idxs = t0 <= obj.t_vect & obj.t_vect < obj.moves_inter(2,law_idx);
             if law_idx == 1
                 x_init = [1,0,0,0;0,0,1,0]*obj.x0;
                 v_init = [0,1,0,0;0,0,0,1]*obj.x0;
             else
                 t_idx_init = find(t_idxs>0,1);
                 x_init = [1,0,0,0;0,0,1,0]*obj.history(:,t_idx_init-1);
                 v_init = [0,1,0,0;0,0,0,1]*obj.history(:,t_idx_init-1);
             end
             syms tf
             x = obj.moves_law(law_idx).law_x(x_init,v_init,t0,tf);
             y = obj.moves_law(law_idx).law_y(x_init,v_init,t0,tf);
             vx = obj.moves_law(law_idx).law_vx(x_init,v_init,t0,tf);
             vy = obj.moves_law(law_idx).law_vy(x_init,v_init,t0,tf);
             
             tf = obj.t_vect(t_idxs);
             x_val = double(subs(x));
             y_val = double(subs(y));
             vx_val = double(subs(vx));
             vy_val = double(subs(vy));
             
             % check dimensions:
             if size(x_val,2)==1
                 x_val = x_val*ones(1,size(tf,2));
             end
             if size(y_val,2)==1
                 y_val = y_val*ones(1,size(tf,2));
             end
             if size(vx_val,2)==1
                 vx_val = vx_val*ones(1,size(tf,2));
             end
             if size(vy_val,2)==1
                 vy_val = vy_val*ones(1,size(tf,2));
             end
             
             obj.history(:,t_idxs)=[x_val;vx_val;y_val;vy_val];
             
          end 
      end
   end
end