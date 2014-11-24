classdef compart     
 methods (Static)
  function J = objective(x)
    Tou = x(1);
    Tu = x(2);
    f=compart;
    inits=[4500 3500 6000];tspan=0:0.1:50;
    [tout,yout]=ode23(@(t,y)f.odesystem(t,y,Tu,Tou),tspan,inits);
     S=yout(end,1);R=yout(end,2);
     N=S+R;%last total value of Tumor 
     U=50*Tou/Tu;%total Angiogenic Inhibitor
     J = N + U;
  end
  function [c,ceq] = constraint(x)
    Tou = x(1);Tu = x(2);
    ceq = [];
    c = Tou-Tu; %nonlinear constraint
  end
  function Dy = odesystem(t,y,Tu,Tou)
        global a c q r gama lamda eta eps v
        at=0;
        for k=0:Tu:50
         at=at + heaviside(t-k)-heaviside(t-(k+Tou)); 
        end 
        u=at;
        S=y(1);R=y(2);K=y(3);
        Dy=[-a*S+(1-v-S/K)*(2-q)*a*S+r*c*R;
        -c*R+(2-r)*c*R*(1-R/K)+(1-v)*q*a*S;
        gama*(S+R)-lamda*K*((S+R)^(2/3))-eta*u*K-eps*v*K];end
 end
end
