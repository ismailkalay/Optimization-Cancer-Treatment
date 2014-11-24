clc;clear;f=hahnfeldt;
%    Tou;Tu   
%    Tou;Tu   
lb=[0;0];
ub=[25;25];
x0=[10,20];
% lb=[0;10];
% ub=[8;15];
% x0=[5,12];
% lb=[0;13];
% ub=[10;25];
% x0=[5,20];
% lb=[0;13];
% ub=[7;25];
% x0=[3,20];
opts = optimset('Algorithm','trust-region-reflective','Display','iter');
% opts = optimset('Algorithm','interior-point','Display','iter');
% opts = optimset('Algorithm','sqp','Display','iter');
% opts = optimset('Algorithm','active-set','Display','iter');
ms = MultiStart('Display','iter');
problem = createOptimProblem('fmincon','x0',x0,'objective',@f.objective,...
'nonlcon',@f.constraint,'lb',lb,'ub',ub,'options',opts);
tic % solution time
% [x,fval,exitflag,output]  = fmincon(@f.objective,x0,[],[],[],[],lb,ub,@f.constraint,opts)
[x,fval,exitflag,output,solutions] = run(ms,problem,3)
toc


