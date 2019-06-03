function a=alphah(v)
global temp
%filename: alphah.m
theta = (v+70)/20;
a=(0.017*exp(0.11*temp))*0.07*exp(-theta);
