function a=alpham(v)
%filename: alpham.m
global temp
theta=(v+45)/10;
if(theta==0)   %check for case that gives 0/0
  a=1.0;  %in that case use L'Hospital's rule
else
  a=(0.017*exp(0.11*temp))*1.0*theta/(1-exp(-theta));
end
