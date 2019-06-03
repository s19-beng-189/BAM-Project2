function a=alphah(v)
global temp
%filename: alphah.m
theta = (v+70)/20;
a=(3^((temp-37)/10))*0.07*exp(-theta);
