clear all;
clc;
x0 = [1 1];
f1 = @(x1,x2) -2*x2+1;
f2 = @(x1,x2) x1-x2;
t0 = 0;
t_f = 150;
h =0.0001;
%% MODIFIED EULER METHOD
[y,t] = modifiedEuler(f1,f2,t0,t_f,x0,h);
k = 0;
t = t0:h:t_f;
for i=1:length((t))
    x1(i)=y(k+1);
    x2(i)=y(k+2);
    k=k+2;
    %disp(x1(i),x2(i))
end
subplot(2,1,1)
plot(t,x1,'LineWidth',2)
legend('x1 by Modified Euler method','x1 by Modified Euler method')
xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
ylabel('x1(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
hold on
subplot(2,1,2)
plot(t,x2,'LineWidth',2)
legend('x2 by Modified Euler method','x2 by Modified Euler method')
xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
ylabel('x2(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
%% MODIFIED EULER METHOD FUNCTION
function [y,t] = modifiedEuler(F1,F2,t0,t_f,y0,h)
y(1) = y0(1);
y(2) = y0(2);
t = t0:h:t_f;
k=0;
for i=1:length((t)-1)
    y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    y1(k+2)=y(k+2)+h*F2(y(k+1),y(k+2));
    y(k+3)=y(k+1)+0.5*h*(F1(y(k+1),y(k+2))+F1(y1(k+1),y1(k+2)));
    y(k+4)=y(k+2)+0.5*h*(F2(y(k+1),y(k+2))+F2(y1(k+1),y1(k+2)));
    k=k+2;
end
end