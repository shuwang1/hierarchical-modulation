function answer=q(x)
syms y;
answer=(1/sqrt(2*pi))*int(exp(-y^2/2),y,x,inf);
answer=eval(answer);
