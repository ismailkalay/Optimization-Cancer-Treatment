clc;clear;
f=compart;%three compartment model
global a c q r gama lamda eta eps v
a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;eta=9.1;
eps=4.7;v=1;
%    Tou;Tu   
% lb = [0;0];%lower-band
% ub = [1;10];%upper-band
% x0 = [0.5,5];% initial guess
% lb=[0;0];
% ub=[4;10];
% x0=[2,6];
lb=[0;0];
ub=[12;25];
x0=[10,15];
% lb=[0;5];
% ub=[3;7];
% x0=[2,6];

opts = psoptimset('Display','iter')

tic % solution time
[xsol,Jval,eflag,outpt] = patternsearch(@f.objective,x0,...
    [],[],[],[],lb,ub,@f.constraint,opts)
toc

