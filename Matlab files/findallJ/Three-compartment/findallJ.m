clear;clc;
% global a c q r gama lamda eta eps 
a=0.02;c=0.2;q=0.9;r=0;gama=5.85;lamda=0.00873;eta=9.1;eps=4.7;
k1=0.5;%Patient resistance to Cancer 
k2=2;%Patient resistance to Chemotherapy 
inits=[4500 3500 6000];tspan=0:1:50;  
counter=0;
tic
for j=1:1:25
    Tu=j;
     for k=1:1:25
     Tou=k;%Tou:open time of U.[days]
     if Tou<=Tu
     [tout,yout]=ode23(@(t,y)compartment3(t,y,Tu,Tou),tspan,inits);
     S=yout(end,1);R=yout(end,2);
     N=S+R;%total last value of Tumor 
     U=50*Tou/Tu;%total Angiogenic Inhibitor
     J(j,k) = k1*N + U*k2;    
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

