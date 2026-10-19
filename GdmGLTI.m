function [A,B,C,D] = GdmGLTI(b,a,Ts,Gdm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
a=[1,a];
b=real(b);
a=real(a);
GLTI=tf(b,a,Ts);

G=GLTI*Gdm;

[A,B,C,D]=ssdata(G);
end