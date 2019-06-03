function b=betam(v)
%filename: betam.m
global temp
theta = (v+70)/18;
b=(0.017*exp(0.11*temp))*4.0*exp(-theta);
