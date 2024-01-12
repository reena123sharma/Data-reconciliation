function K_e = K_h (x,A_r,h)
  K_e = inv(det(A))*(h^(-1))*K_x((inv(A_r))*x/h)
endfunction