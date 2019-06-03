function b=betan(v)
%filename: betan.m
global temp
theta = (v+70)/80;
b=(3^((temp-37)/10))*0.125*exp(-theta);
