clc;
clear all;
% T = time 
% u0 = initial conditions
% N = number of intervals
% K = number of iterations
% U = error (initially 0)
% MG = no. of coarse intervals
% MF = no. of fine intervals
T1=50;
u01=[-10 -10 0.1 0.1 0.1 0.1 0.1 0.1];
N1=1000;
K1=1;
MG1=1;
MF1=10;
U1 = zeros(N1,N1+1,8);
parareal(T1, u01, N1, K1, U1, MG1, MF1);

function parareal(T, u0, N, K, U, MG, MF)
    dT = T / N;
    TT = (0:N) * dT;
    Go = zeros(N+1, length(u0));
    Gn = zeros(N+1, length(u0));
    Fn = zeros(N+1, length(u0));
    
    U(1, 1, :) = u0;
    
    for n = 1:N
        Go(n + 1, :) = Coarse(TT(n), TT(n + 1), squeeze(U(1, n, :)), MG);
        U(1, n + 1, :) = Go(n + 1, :);
        n=n+1;
    end
    A = 8;
    B = 1;
    
    for k = 1:K
        for n = 1:A
            for m = 1:B
                % OpenMP-like parallel processing using parfor
                Fine(TT(m + B*(n - 1)), TT(m + B*n), squeeze(U(k, m + B*(n - 1), :)), MF);
            end
        end
    end
    U(k + 1, 1, :) = u0;
    
    for n = 1:N
        Gn(n + 1, :) = Coarse(TT(n), TT(n + 1), squeeze(U(k + 1, n, :)), MG);
        for i = 1:length(u0)
            U(k + 1, n + 1, i) = Fn(n + 1, i) + Gn(n + 1, i) - Go(n + 1, i);
        end
    end

     k = 0;
    h = T/N*10;
    t = 0:h:T;
    for i=1:length((t))
         x1(i)=Go(k+1);
         x2(i)=Go(k+2);
         x3(i)=Go(k+3);
         x4(i)=Go(k+4);
         x5(i)=Go(k+5);
         x6(i)=Go(k+6);
         x7(i)=Go(k+7);
         x8(i)=Go(k+8);
         k=k+8;
    end
    plot(t,x1)
    subplot(2,1,1)
    plot(t,x1,'LineWidth',2)
    legend('x1 by Parareal method','x1 by Parareal method')
    xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    ylabel('x1(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    hold on
    subplot(2,1,2)
    plot(t,x2,'LineWidth',2)
    legend('x2 by Parareal method','x2 by Parareal method')
    xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    ylabel('x2(t)','FontSize',16,'FontName','Arial','FontWeight','bold')

    for n = 1:N
        for i = 1:length(u0)
            Go(n, i) = Gn(n, i);
            %disp(Go)
        end
    end

    % k = 0;
    % h = T/N*10;
    % t = 0:h:T;
    % for i=1:length((t))
    %      x1(i)=Go(k+1);
    %      x2(i)=Go(k+2);
    %      x3(i)=Go(k+3);
    %      x4(i)=Go(k+4);
    %      x5(i)=Go(k+5);
    %      x6(i)=Go(k+6);
    %      x7(i)=Go(k+7);
    %      x8(i)=Go(k+8);
    %      k=k+8;
    % end
    % plot(t,x1)
    % subplot(2,1,1)
    % plot(t,x1,'LineWidth',2)
    % legend('x1 by Parareal method','x1 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x1(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % hold on
    % subplot(2,1,2)
    % plot(t,x2,'LineWidth',2)
    % legend('x2 by Parareal method','x2 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x2(t)','FontSize',16,'FontName','Arial','FontWeight','bold')

    % subplot(8,1,3)
    % plot(t,x3,'LineWidth',2)
    % legend('x3 by Parareal method','x3 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x3(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,4)
    % plot(t,x4,'LineWidth',2)
    % legend('x4 by Parareal method','x4 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x4(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,5)
    % plot(t,x5,'LineWidth',2)
    % legend('x5 by Parareal method','x5 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x5(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % hold on
    % subplot(8,1,6)
    % plot(t,x6,'LineWidth',2)
    % legend('x6 by Parareal method','x6 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x6(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,7)
    % plot(t,x7,'LineWidth',2)
    % legend('x7 by Parareal method','x7 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x7(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    % subplot(8,1,8)
    % plot(t,x8,'LineWidth',2)
    % legend('x8 by Parareal method','x8 by Parareal method')
    % xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
    % ylabel('x8(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
    
end

function Go = Coarse(tn, tn1, un, MG) %function [y,t] = euler(t0,t_f,y0,h);
    ode1=@(x1,x2,x3,x4,x5,x6,x7,x8) x1;
    ode2=@(x1,x2,x3,x4,x5,x6,x7,x8) x2;
    ode3=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
    ode4=@(x1,x2,x3,x4,x5,x6,x7,x8) x4;
    ode5=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
    ode6=@(x1,x2,x3,x4,x5,x6,x7,x8) x6;
    ode7=@(x1,x2,x3,x4,x5,x6,x7,x8) x7;
    ode8=@(x1,x2,x3,x4,x5,x6,x7,x8) x8;
    % Implement your G function here
    % Example: Go = YourGFunction(tn, tn1, un, MG);
    h = (tn1-tn)/MG;
    Go(1)=un(1); %y(1) = y0(1);
    Go(2)=un(2); %y(2) = y0(2);
    Go(3)=un(3); 
    Go(4)=un(4);
    Go(5)=un(5); 
    Go(6)=un(6);
    Go(7)=un(7); 
    Go(8)=un(8);
    t=tn:h:tn1;  %t = t0:h:t_f;
    l=0; %k=0;
    for i=1:length((t)-1) %for i=1:length((t)-1)
        Go(l+9)=Go(l+1)+h*ode1(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));%y(k+3)=y(k+1)+h*F1(y(k+1),y(k+2));
        Go(l+10)=Go(l+2)+h*ode2(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+11)=Go(l+3)+h*ode3(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+12)=Go(l+4)+h*ode4(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+13)=Go(l+5)+h*ode5(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+14)=Go(l+6)+h*ode6(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+15)=Go(l+7)+h*ode7(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        Go(l+16)=Go(l+8)+h*ode8(Go(l+1),Go(l+2),Go(l+3),Go(l+4),Go(l+5),Go(l+6),Go(l+7),Go(l+8));
        l=l+8; %k=k+8;
    end
    Go = [Go(2*MG + 1), Go(2*MG + 2), Go(2*MG + 3), Go(2*MG + 4), Go(2*MG + 5), Go(2*MG + 6), Go(2*MG + 7), Go(2*MG + 8)];
end

 function [Fn,t]=Fine(tn, tn1, un, MF) %function [y,t] = modifiedEuler(t0,t_f,y0,h);
    % Implement your F function here
    % Example: Fn = YourFFunction(tn, tn1, un, MF);
    ode1=@(x1,x2,x3,x4,x5,x6,x7,x8) x1;
    ode2=@(x1,x2,x3,x4,x5,x6,x7,x8) x2;
    ode3=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
    ode4=@(x1,x2,x3,x4,x5,x6,x7,x8) x4;
    ode5=@(x1,x2,x3,x4,x5,x6,x7,x8) 0;
    ode6=@(x1,x2,x3,x4,x5,x6,x7,x8) x6;
    ode7=@(x1,x2,x3,x4,x5,x6,x7,x8) x7;
    ode8=@(x1,x2,x3,x4,x5,x6,x7,x8) x8;

    h = (tn1-tn)/MF;

    Fn(1)=un(1); %y(1) = y0(1);
    Fn(2)=un(2); %y(2) = y0(2);
    Fn(3)=un(3); 
    Fn(4)=un(4);
    Fn(5)=un(5); 
    Fn(6)=un(6);
    Fn(7)=un(7); 
    Fn(8)=un(8);
    t = tn:h:tn1; %t = t0:h:t_f;
    l=0; %k=0;
    for i=1:length((t)-1) %for i=1:length((t)-1)
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
    Fn = [Fn(2*MF+1),Fn(2*MF+2),Fn(2*MF+3),Fn(2*MF+4),Fn(2*MF+5),Fn(2*MF+6),Fn(2*MF+7),Fn(2*MF+8) ];
end
