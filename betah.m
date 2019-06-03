function b=betah(v)
%filename: betah.m
global t
theta = (v+40)/10;
b=1.0/(1+exp(-theta));