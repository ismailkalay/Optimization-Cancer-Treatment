classdef continuous    
 methods (Static)
    function Dy = hahnfeldtmodel(t,y)
        % constants hahnfeldt in A.Swierniak 2013
        beta=0.192;gama=5.85;lamda=0.00873;
        mu=0;eta=0.15;eps=0.26;psi=0.34;u = heaviside(t);v=heaviside(t);
        N=y(1);K=y(2);
        Dy=[-beta*N*log(N/K)-psi*v*N;
            gama*N-lamda*K*(N^(2/3))-mu*K-eta*u*K-eps*v*K];end
    function Dy = hpruningmodel(t,y)
        % constants in A.Swierniak 2013
        beta=0.192;gama=5.85;lamda=0.00873;
        mu=0;eta=0.15;eps=0.26;psi=0.34;u=20;v=2;N=y(1);K=y(2);
%         psi=(0.3/(1+(((K/N)-2)/0.35)^2));
        Dy=[-beta*N*log(N/K)-(0.3/(1+(((K/N)-2)/0.35)^2))*v*N;
            gama*K-lamda*K*(N^(2/3))-mu*K-eta*u*K-eps*v*K];end
    function Dy = benzekrymodel(t,y)
        % constants Benzekry in A.Swierniak 2013
        beta=0.192;gama=5.85;lamda=0.00873;
        eps=7.56e-3;tou=7.5e-3;eta=6.85e-7;psi=1.37e-5;
        N=y(1);M=y(2);I=y(3);Q=(M/(M+I));
        u=525/2*heaviside(t);
        v=(212/2)/7*heaviside(t);
%          at=0;Tu=7;Gu=1;
%         for k=0:Tu:50
%          at=at + heaviside(t-k)-heaviside(t-(k+Gu)); 
%         end 
%         v=(212/2)*at;
        Dy=[-beta*N*log(N/M)-psi*v*N*Q*M;
            eps*I-tou*M;
            -eps*I+gama*N-lamda*I*(N^(2/3))-eta*u*I*Q*M];end
    function Dy = compartmentmodel(t,y)
        % constants in A.SwierSiak 2013
        a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;eta=9.1;eps=4.7;
        S=y(1);R=y(2);K=y(3);v=1;u=1;
        Dy=[-a*S+(1-v-S/K)*(2-q)*a*S+r*c*R;
            -c*R+(2-r)*c*R*(1-R/K)+(1-v)*q*a*S;
            gama*(S+R)-lamda*K*((S+R)^(2/3))-eta*u*K-eps*v*K];end
    
    function hahnfeldt
    inits=[8000 10000];tspan=0:0.1:50;f=continuous;
   
    [tout,yout]=ode23(@f.hahnfeldtmodel,tspan,inits);
    f.PlotHahnfeldt(tout,yout)
    end
    function hpruning
    inits=[8000 10000];tspan=[0 50];f=continuous;  
    [tout,yout]=ode23(@f.hpruningmodel,tspan,inits);
    f.PlotHPruning(tout,yout)
    end
    function benzekry
    inits=[7900 2000 6000];tspan=[0 50];f=continuous;  
    [tout,yout]=ode23(@f.benzekrymodel,tspan,inits);
    f.PlotBenzekry(tout,yout)
    end
    function out = compartment
    inits=[4500 3500 6000];tspan=[0 50];f=continuous;
    [tout,yout]=ode23(@f.compartmentmodel,tspan,inits);
    f.PlotCompartment(tout,yout);
    S=yout(end,1);R=yout(end,2);K=yout(end,3);
    N=S+R;%last value of Tumor 
    out =[N,S,R,K];   
    end

    function PlotHahnfeldt(T,Y)
%     figure
%     subplot(2,2,1)
    plot(T,Y(:,1),'-',T,Y(:,2),'-');
    title('Hahnfeldt Model');xlabel('Time[Days]');ylabel('Volume[mm^3]');
    legend('N:Cancer Volume','K:Endothelial Volume');end
    function PlotHPruning(T,Y)
%     subplot(2,2,2)
    plot(T,Y(:,1),'-',T,Y(:,2),'-');
    title('Hahnfeldt Pruning Model');xlabel('Time[Days]');ylabel('Volume[mm^3]');
    legend('N:Cancer Volume','K:Endothelial Volume');end
    function PlotBenzekry(T,Y)
		plot(T,Y(:,1),'-',T,Y(:,2),'-',T,Y(:,3),'-');
		title('benzekry Model');
		xlabel('Time[Days]');ylabel('Volume[mm^3]');
		legend('N:Cancer Volume','M:Mature Vessels',...
		'I:Immature Vessels');end  
    function PlotCompartment(T,Y)
%     subplot(2,2,4)
    plot(T,Y(:,1),'-',T,Y(:,2),'-',T,Y(:,3),'-');
    title('Tree-compartment model');xlabel('Time[Days]');ylabel('Volume[mm^3]');
    legend('S:Sensitive Cancer','R:Resistant Cancer','K:Endothelial Volume');end    
 end
end
 
