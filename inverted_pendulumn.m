function [] = simsprpen()
clc;clf;clear;
p = setpar();
tspan = 10;
dt = 0.1;
N = (tspan-0)/dt;
L = 1;
%%...........runge-kutta...............%%
[t1,z1] = ode23(@springpend,[0 tspan],[1;0;pi/4;0],[],p);
subplot(2,1,1)
plot(t1,z1(:,1),t1,z1(:,2),t1,z1(:,3),t1,z1(:,4));
hold on;
%%...........Dormand-prince............%%
[t2,z2] = ode45(@springpend,[0 tspan],[1;0;pi/4;pi/4],[],p);
plot(t2,z2(:,1),'--',t2,z2(:,2),'--',t2,z2(:,3),'--',t2,z2(:,4),'--');
legend('dx','x','d\phi','\phi');
%%...........Analysis of differences.....%%
diffy = z1(:,1)-interp1(t2,z2(:,1),t1);
%%...........Comparison between methods...%%
title('Runge-kutta vs Dormand-Prince');
subplot(2,2,3);
plot(t1,z1(:,1),t2,z2(:,1),t1,diffy);
xlabel('t[sec]');ylabel('x[m]');
title('For different initial values');
legend('x1','x2','x1-x2');
subplot(2,2,4);
plot(z1(:,1),z1(:,2),'r',z2(:,1),z2(:,2),'b');
xlabel('v[m/sec]');
ylabel('x[m]');
title('Runge-kutta vs Dormand-Prince');
grid;
%%...........XY plot of Spring Pendulumn..........%%
del = zeros(N, 1);
theta = zeros(N, 1);
x = zeros(N, 1);
y = zeros(N, 1);
for i = 1:N
    theta(i) = z1(i,3);
    del(i) = z1(i,1);
    x(i) = L*(+del(i)*sin(theta(i)));
    y(i) = L*(-(del(i))*cos(theta(i)));
    pause(0.01);
    figure(2);
    axis([-1,1,-2,0]);
    plot(x, y,'r');
    hold on;
    title('XY Plot');
    xlabel('x position');
    ylabel('y position');
end
%%............Parameters..................%%
function [p]= setpar()
p.g =9.8;
p.m = 0.1;
p.C = 3;
p.Lo = 1;
%%...........Defining System...............%%
function yp = springpend(t,z,p)
LL = z(1);
vLL=z(2);
A=z(3);
vA=z(4);
dLL = vLL;
dvLL = (p.g*p.m*cos(A)-p.C*LL+p.C*p.Lo+p.m*LL*vA^2)/p.m;
dA = vA;
dvA = -(p.g*sin(A)+(2)*(vLL*vA))/LL;
yp = [dLL;dvLL;dA;dvA];

