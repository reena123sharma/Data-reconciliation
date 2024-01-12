function K_e = K_x (x)
  if norm(x) < 1
    K_e = (2/pi)*(1-norm(x)^2);
   else
    K_e = 0;
  end
  end