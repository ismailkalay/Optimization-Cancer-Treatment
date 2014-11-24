gs = GlobalSearch;
sixmin = @(x)(4*x(1)^2 - 2.1*x(1)^4 + x(1)^6/3 ...
    + x(1)*x(2) - 4*x(2)^2 + 4*x(2)^4);
problem = createOptimProblem('fmincon','x0',[-1,2],...
    'objective',sixmin,'lb',[-3,-3],'ub',[3,3]);
[xmin,fmin,flag,outpt,allmins] = run(gs,problem)

options = optimset('Algorithm','sqp','Display','iter');
% options = optimset('Algorithm','trust-region-reflective','Display','iter');
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',@f.objective,'lb',lb,'ub',ub,'options',options);
tic % solution time
% [x,fval,exitflag,output]  = fmincon(@f.objective,x0,[],[],[],[],lb,ub,@f.constraint,options)
[x,fval] = run(GlobalSearch,problem)
toc