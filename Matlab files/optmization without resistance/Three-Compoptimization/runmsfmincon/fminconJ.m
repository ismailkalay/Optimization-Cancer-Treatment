clc;clear;f=compart;
%three compartment model
global a c q r gama lamda eta eps v
a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;eta=9.1;eps=4.7;v=1;
%    Tou;Tu   
% lb = [0;0];%lower-band
% ub = [1;10];%upper-band
% x0 = [0.5,5];% initial guess
% lb=[0;0];
% ub=[4;10];
% x0=[2,6];
lb=[0;0];
ub=[25;25];
x0=[5,15];
% lb=[0;5];
% ub=[3;7];
% x0=[2,6];
% opts = optimset('GradObj','off');
% opts = optimset('Algorithm','trust-region-reflective');
% opts = optimset('Algorithm','active-set','Display','iter');
% opts = optimset('Algorithm','interior-point','Display','iter');
% opts = optimset('Algorithm','sqp','Display','iter');
% ms = MultiStart('Display','iter');
% problem = createOptimProblem('fmincon','x0',x0,'objective',@f.objective,...
% 'nonlcon',@f.constraint,'lb',lb,'ub',ub,'options',opts);
tic % solution time
[x,fval,exitflag,output]  = fmincon(@f.objective,x0,[],[],[],[],lb,ub,@f.constraint,opts)
% [x,fval,exitflag,output,solutions] = run(ms,problem,1)
toc


