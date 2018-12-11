function [new_pt] = new_dist_angle(pt, a, l) 
 
    new_pt = pt + [cos(a)*l; sin(a)*l; 0];
end
 

