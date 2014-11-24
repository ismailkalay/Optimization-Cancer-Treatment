function Dy = compartment3(t,y,Tu,Tou)
        % constants in A.SwierSiak 2013    
         at=0;
        for k=0:Tu:50
         at=at + heaviside(t-k)-heaviside(t-(k+Tou)); 
        end 
        beta=0.192;gama=5.85;lamda=0.00873;
        eps=7.56e-3;tou=7.5e-3;eta=6.85e-7;psi=1.37e-5;
        N=y(1);M=y(2);I=y(3);Q=(M/(M+I));
        u=(525/2)*at;
        v=(212/2)/7*heaviside(t);
        Dy=[-beta*N*log(N/M)-psi*v*N*Q*M;
            eps*I-tou*M;
            -eps*I+gama*N-lamda*I*(N^(2/3))-eta*u*I*Q*M];end
        
    
    
    

       