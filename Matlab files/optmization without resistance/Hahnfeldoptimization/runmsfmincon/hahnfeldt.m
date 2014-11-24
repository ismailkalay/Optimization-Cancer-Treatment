classdef hahnfeldt     
 methods (Static)
  function J = objective(x)
    Tou = x(1);
    Tu = x(2);
    f=hahnfeldt;k1=0.5;k2=2;
    inits=[8000 10000];tspan=0:0.1:50;
    [tout,yout]=ode23(@(t,y)f.odesystem(t,y,Tu,Tou),tspan,inits);
     N=yout(end,1); %total last value of Tumor 
     U=50*Tou/Tu; %total Angiogenic Inhibitor
     J = k1*N + k2*U; 
  end
  function [c,ceq] = constraint(x)
    Tou = x(1);Tu = x(2);
    ceq = [];
    c = Tou-Tu; %nonlinear constraint
  end
  function Dy = odesystem(t,y,Tu,Tou)
       at=0;
        for k=0:Tu:50
         at=at + heaviside(t-k)-heaviside(t-(k+Tou)); 
        end 
        beta=0.192;gama=5.85;lamda=0.00873;
        mu=0;eta=0.15;eps=0.26;psi=0.34;u = at;v=heaviside(t);
        N=y(1);K=y(2);
        Dy=[-beta*N*log(N/K)-psi*v*N;
            gama*N-lamda*K*(N^(2/3))-mu*K-eta*u*K-eps*v*K];end
 end
end