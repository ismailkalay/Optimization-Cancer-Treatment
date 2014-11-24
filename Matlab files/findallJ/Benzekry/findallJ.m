clear;clc;
k1=0.5;%Patient resistance to Cancer 1<k1<2
k2=2;%Patient resistance to Chemotherapy 1<k2<2
inits=[7900 2000 6000];tspan=0:0.1:50;counter=0;
tic
for j=1:1:25
    Tu=j;
     for k=1:1:25
     Tou=k;%Tou:open time of U.[days]
     if Tou<=Tu
     [tout,yout]=ode23(@(t,y)benzekry(t,y,Tu,Tou),tspan,inits);
      N=yout(end,1); %total last value of Tumor 
     U=50*Tou/Tu; %total Angiogenic Inhibitor
     J(j,k) = k1*N + k2*U;    
     else
         counter = counter + 1;
     end
     end
end
toc
figure(1)
mesh(J);
colormap(hsv)  
xlabel('Tou')
ylabel('Tu')
zlabel('J')

% figure(2)
% surfl(J);
% shading interp;
% colormap(pink);
% xlabel('Tou')
% ylabel('Tu')
% zlabel('J')

% tol = 1.e-5;
% J(abs(imag(J))<tol) = real(J(abs(imag(J))<tol));

