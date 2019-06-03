function a=alphan(v)
%filename: alphan.m
global temp
theta=(v+60)/10;
if(theta==0)   %check for case that gives 0/0
  a=0.1;  %in that case use L'Hospital's rule
else
  a=((3^((temp-37)/10))*0.1*theta)/(1-exp(-theta));
end
