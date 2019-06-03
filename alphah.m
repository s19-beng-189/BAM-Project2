function a=alphah(v)
%filename: alphah.m
global t 
theta = (v+70)/20;
a=0.07*exp(-theta);