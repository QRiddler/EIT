clc;
clear all;
global ODE1 ODE2
ODE1 = @(x1,x2) -2*x2+1;
ODE2=@(x1,x2) x1-x2;
% T = time 
% u0 = initial conditions
% N = number of intervals
% K = number of iterations
% U = error (initially 0)
% MG = no. of coarse intervals
% MF = no. of fine intervals

T1=150;
u01=[1 0];
N1=10000;
K1=90;
MG1=10;
MF1=50;
U1 = zeros(N1,N1+1,2);
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
    A = 100;
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
         k=k+2;
    end
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
    %      k=k+2;
    % end
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
end

function [Go,t] = Coarse(tn, tn1, un, MG) %function [y,t] = euler(t0,t_f,y0,h);
    ODE1=@(x1,x2) -2*x2+1;
    ODE2=@(x1,x2) x1-x2;
    % Implement your G function here
    % Example: Go = YourGFunction(tn, tn1, un, MG);
    h = (tn1-tn)/MG;
    Go(1)=un(1); %y(1) = y0(1);
    Go(2)=un(2); %y(2) = y0(2);
    t=tn:h:tn1;  %t = t0:h:t_f;
    l=0; %k=0;
    for i=1:length((t)-1) %for i=1:length((t)-1)
        Go(l+3)=Go(l+1)+h*ODE1(Go(l+1),Go(l+2));  %y(k+3)=y(k+1)+h*F1(y(k+1),y(k+2));
        Go(l+4)=Go(l+2)+h*ODE2(Go(l+1),Go(l+2)); %y(k+4)=y(k+2)+h*F2(y(k+1),y(k+2));
        l=l+2; %k=k+2;
    end
    Go = [Go(2*MG+1), Go(2*MG + 2)];
 end

 function [Fn,t]=Fine(tn, tn1, un, MF) %function [y,t] = modifiedEuler(t0,t_f,y0,h);
    % Implement your F function here
    % Example: Fn = YourFFunction(tn, tn1, un, MF);
    ODE1=@(x1,x2) -2*x2+1;
    ODE2=@(x1,x2) x1-x2;
    h = (tn1-tn)/MF;
    Fn(1)=un(1); %y(1) = y0(1);
    Fn(2)=un(2); %y(2) = y0(2);
    t = tn:h:tn1; %t = t0:h:t_f;
    l=0; %k=0;
    for i=1:length((t)-1) %for i=1:length((t)-1)
    var(l+1)=Fn(l+1)+h*ODE1(Fn(l+1),Fn(l+2)); %y1(k+1)=y(k+1)+h*F1(y(k+1),y(k+2));
    var(l+2)=Fn(l+2)+h*ODE2(Fn(l+1),Fn(l+2)); %y1(k+2)=y(k+2)+h*F2(y(k+1),y(k+2));
    Fn(l+3)=Fn(l+1)+0.5*h*(ODE1(Fn(l+1),Fn(l+2))+ODE1(var(l+1),var(l+2))); %y(k+3)=y(k+1)+0.5*h*(F1(y(k+1),y(k+2))+F1(y1(k+1),y1(k+2)));
    Fn(l+4)=Fn(l+2)+0.5*h*(ODE2(Fn(l+1),Fn(l+2))+ODE2(var(l+1),var(l+2))); %y(k+4)=y(k+2)+0.5*h*(F2(y(k+1),y(k+2))+F2(y1(k+1),y1(k+2)));
    l=l+2; %k=k+2;
    end
    Fn = [Fn(2*MF + 1),Fn(2*MF+2)];
end
