clc;clear;
f=hahnfeldt;%three compartment model
%    Tou;Tu   
lb=[0;0];
ub=[25;25];
x0=[10,15];
% lb=[0;10];
% ub=[8;15];
% x0=[5,12];
% lb=[0;13];
% ub=[10;25];
% x0=[5,20];
% lb=[0;13];
% ub=[7;25];
% x0=[3,20];

opts = psoptimset('Display','iter');
tic % solution time
[xsol,Jval,eflag,outpt] = patternsearch(@f.objective,x0,...
    [],[],[],[],lb,ub,@f.constraint,opts)
toc

