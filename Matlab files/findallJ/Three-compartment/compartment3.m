function Dy = compartment3(t,y,Tu,Tou)
        % constants in A.SwierSiak 2013    
         at=0;
        for k=0:Tu:50
         at=at + heaviside(t-k)-heaviside(t-(k+Tou)); 
        end 
        global a c q r gama lamda eta eps 
%         a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;
% 		eta=9.1;eps=4.7;
		S=y(1);R=y(2);K=y(3);v=heaviside(t);u=at;
		Dy=[-a*S+(1-v-S/K)*(2-q)*a*S+r*c*R;
            -c*R+(2-r)*c*R*(1-R/K)+(1-v)*q*a*S;
            gama*(S+R)-lamda*K*((S+R)^(2/3))-eta*u*K-eps*v*K];end

