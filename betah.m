function b=betah(v)
%filename: betah.m
global temp
theta = (v+40)/10;
b=((0.017*exp(0.11*temp))*1.0)/(1+exp(-theta));
