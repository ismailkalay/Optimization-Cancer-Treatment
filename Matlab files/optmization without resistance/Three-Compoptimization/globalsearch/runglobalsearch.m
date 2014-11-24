clear;clc;f=compart;
%three compartment model
 global a c q r gama lamda eta eps v
a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;eta=9.1;
eps=4.7;v=1;
%    Tou;Tu   
% lb = [0;0];%lower-band
% ub = [1;10];%upper-band
% x0 = [0.5,5];% initial guess
% lb=[0;0];
% ub=[4;10];
% % x0=[2,6];
% lb=[0;0];
% ub=[25;25];
% x0=[10,15];
lb=[0;5];
ub=[3;7];
x0=[2,6];
% opts = optimset('GradObj','off');
% opts = optimset('Algorithm','active-set');
% opts = optimset('Algorithm','sqp','Display','iter');
% opts = optimset('Display','iter');
% opts = optimset('Algorithm','trust-region-reflective','Display','iter');
gs = GlobalSearch('Display','iter');
opts = optimset('Algorithm','active-set','Display','iter');
problem = createOptimProblem('fmincon','x0',x0,'objective',@f.objective,...
'nonlcon',@f.constraint,'lb',lb,'ub',ub,'options',opts);
[xmin,fmin,flag,outpt,allmins] = run(gs,problem);



