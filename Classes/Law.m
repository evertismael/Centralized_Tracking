classdef Law
   properties
      law_x
      law_y
      law_vx
      law_vy
   end
   methods
      function obj = Law(law_x,law_y,law_vx, law_vy)
          obj.law_x = law_x;
          obj.law_y = law_y;
          obj.law_vx = law_vx;
          obj.law_vy = law_vy;
      end
   end
end