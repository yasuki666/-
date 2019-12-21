clear; 
%参数设置 
Re=10; %雷诺数取10,100，500，1000 
L=1; %空穴几何尺寸 
n=100; 
dh=L/n;%delta h 
dt=1e-4; %时间步长 
psi=zeros(n+1,n+1); 
xi=zeros(n+1,n+1); 
rho=1; 
for k=1:1000000 
err=0; 
%边界条件 
for i=2:n 
xi(i,1)=-2*(psi(i,2)-psi(i,1))/dh^2; 
xi(i,n+1)=-2*(psi(i,n)-psi(i,n+1))/dh^2; 
end 
for j=2:n 
xi(1,j)=-2*(psi(2,j)-psi(1,j)+dh)/dh^2; 
xi(n+1,j)=-2*(psi(n,j)-psi(n+1,j))/dh^2; 
end 
%控制方程 
for i=2:n 
for j=2:n 
u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*dh); 
v(i,j)=-((psi(i+1,j)-psi(i-1,j))/(2*dh)); 
err1=(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+xi(i,j)*dh^2)/4-psi(i,j); 
psi(i,j)=psi(i,j)+rho*err1; 
err2=dt*(-dh/2*(u(i,j)*(xi(i+1,j)-xi(i-1,j))+v(i,j)*(xi(i,j+1)-xi(i,j-1)))+(xi(i+1,j)+xi(i-1,j)+xi(i,j+1)+xi(i,j-1)-4*xi(i,j))/Re)/dh^2; 
xi(i,j)=xi(i,j)+rho*err2; 
temp=max(abs(err1),abs(err2)); 
if err<temp 
err=temp; 
end 
end 
end 
if (mod(k,1000)==0) %每千步显示结果 
k 
err 
contour(psi,100);%contour求迹线 
pause(0.5) 
end 
if err<1e-6 
break; 
end 
end 
k 
err 
rho 
dt 
contour(psi,100); 
