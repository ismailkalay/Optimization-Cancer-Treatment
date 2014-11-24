function Dy = compartment3(t,y,Tu,Tou)
        % constants in A.SwierSiak 2013    
         at=0;
        for k=0:Tu:50
         at=at + heaviside(t-k)-heaviside(t-(k+Tou)); 
        end 
%         global beta psi mu gama lamda eta eps 
% constants hahnfeldt in A.Swierniak 2013
        beta=0.192;gama=5.85;lamda=0.00873;
        mu=0;eta=0.15;eps=0.26;psi=0.34;u = at;v=heaviside(t);
        N=y(1);K=y(2);
        Dy=[-beta*N*log(N/K)-psi*v*N;
            gama*N-lamda*K*(N^(2/3))-mu*K-eta*u*K-eps*v*K];end
        