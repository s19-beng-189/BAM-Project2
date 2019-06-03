function b=betam(v)
%filename: betam.m
global temp
theta = (v+70)/18;
b=(3^((temp-37)/10))*4.0*exp(-theta);
