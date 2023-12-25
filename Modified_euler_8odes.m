clear all;
clc;
x0 = [-10 -10 0.1 0.1 0.1 0.1 0.1 0.1];
 ode1=@(x1,x2,x3,x4,x5,x6,x7,x8) x1;
 ode2=@(x1,x2,x3,x4,x5,x6,x7,x8) x2;
 ode3=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
 ode4=@(x1,x2,x3,x4,x5,x6,x7,x8) x4;
 ode5=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
 ode6=@(x1,x2,x3,x4,x5,x6,x7,x8) x6;
 ode7=@(x1,x2,x3,x4,x5,x6,x7,x8) x7;
 ode8=@(x1,x2,x3,x4,x5,x6,x7,x8) x8;
t0 = 0;
t_f = 50;
h =0.0001;
%% MODIFIED EULER METHOD
[y,t] = modifiedEuler(ode1,ode2,ode3,ode4,ode5,ode6,ode7,ode8,t0,t_f,x0,h);
k = 0;
t = t0:h:t_f;
for i=1:length((t))
    x1(i)=y(k+1);
    x2(i)=y(k+2);
    x3(i)=y(k+3);
    x4(i)=y(k+4);
    x5(i)=y(k+5);
    x6(i)=y(k+6);
    x7(i)=y(k+7);
    x8(i)=y(k+8);
    k=k+8;
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
    % subplot(8,1,1)
    % plot(t,x1,'LineWidth',2)
    % legend('x1 by Modified Euler method','x1 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x1(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % hold on
    % subplot(8,1,2)
    % plot(t,x2,'LineWidth',2)
    % legend('x2 by Modified Euler method','x2 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x2(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,3)
    % plot(t,x3,'LineWidth',2)
    % legend('x3 by Modified Euler method','x3 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x3(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,4)
    % plot(t,x4,'LineWidth',2)
    % legend('x4 by Modified Euler method','x4 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x4(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,5)
    % plot(t,x5,'LineWidth',2)
    % legend('x5 by Modified Euler method','x5 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x5(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % hold on
    % subplot(8,1,6)
    % plot(t,x6,'LineWidth',2)
    % legend('x6 by Modified Euler method','x6 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x6(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,7)
    % plot(t,x7,'LineWidth',2)
    % legend('x7 by Modified Euler method','x7 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x7(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,8)
    % plot(t,x8,'LineWidth',2)
    % legend('x8 by Modified Euler method','x8 by Modified Euler method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x8(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
%% MODIFIED EULER METHOD FUNCTION
function [Fn,t] = modifiedEuler(ode1,ode2,ode3,ode4,ode5,ode6,ode7,ode8,t0,t_f,y0,h)
    Fn(1)=y0(1); %y(1) = y0(1);
    Fn(2)=y0(2); %y(2) = y0(2);
    Fn(3)=y0(3); 
    Fn(4)=y0(4);
    Fn(5)=y0(5); 
    Fn(6)=y0(6);
    Fn(7)=y0(7); 
    Fn(8)=y0(8);
    t = t0:h:t_f;
    l=0;
for i=1:length((t)-1)
    var(l+1)=Fn(l+1)+h*ode1(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8)); %y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    var(l+2)=Fn(l+2)+h*ode2(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8)); %y1(k+2)=y(k+2)+h*F2(y(k+1),y(k+2));
    var(l+3)=Fn(l+3)+h*ode3(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8)); %y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    var(l+4)=Fn(l+4)+h*ode4(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8));
    var(l+5)=Fn(l+5)+h*ode5(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8)); %y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    var(l+6)=Fn(l+6)+h*ode6(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8));
    var(l+7)=Fn(l+7)+h*ode7(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8)); %y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    var(l+8)=Fn(l+8)+h*ode8(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8));
    Fn(l+9)=Fn(l+1)+0.5*h*(ode1(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode1(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8))); %y(k+3)=y(k+1)+0.5*h*(F1(y(k+1),y(k+2))+F1(y1(k+1),y1(k+2)));
    Fn(l+10)=Fn(l+2)+0.5*h*(ode2(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode2(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+11)=Fn(l+3)+0.5*h*(ode3(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode3(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+12)=Fn(l+4)+0.5*h*(ode4(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode4(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+13)=Fn(l+5)+0.5*h*(ode5(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode5(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+14)=Fn(l+6)+0.5*h*(ode6(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode6(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+15)=Fn(l+7)+0.5*h*(ode7(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode7(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    Fn(l+16)=Fn(l+8)+0.5*h*(ode8(Fn(l+1),Fn(l+2),Fn(l+3),Fn(l+4),Fn(l+5),Fn(l+6),Fn(l+7),Fn(l+8))+ode8(var(l+1),var(l+2),var(l+3),var(l+4),var(l+5),var(l+6),var(l+7),var(l+8)));
    l=l+8; %k=k+8;
end
end