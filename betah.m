function b=betah(v)
%filename: betah.m
global temp
theta = (v+40)/10;
b=((3^((temp-37)/10))*1.0)/(1+exp(-theta));
