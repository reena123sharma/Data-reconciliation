function y = pdf (x,Np,w,x_k,A,h)
  sum = 0;
  for i=1:Np 
    sum = sum +  w(i)*inv(det(A))*(h^(-1))*K_x((inv(A))*(x - x_k(i))/h);
  end   
  y = sum;
  end