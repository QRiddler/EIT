clc;
clear all;
global Xe omega_s v_dr  v_qr  bk
global v Vb omega_elB omega_el Kmrr R2 R1 Ls_prime Tr Pt k c Ht Hg Xm Rs Rr Lm Lss Lrr

%% measurements are taken considering the backlash amount equal to 0.2
bk=0.1;
% Xe=.06;
% Xm = 4;
%% Xm and Xe data are adjusted to find the transient free Qe and Pe
Xe=.01;
Xm = 6;
Rs = Xm/800;
Rr = 1.1*Rs;
Lm = Xm;
Lss = 1.01*Lm;
Lrr = 1.005*Lss;
omega_s = 1;
omega_elB = 314;
k = 0.3; 
c = 0.01;
Ht = 4; 
Hg = 0.1*Ht;
omega_el = omega_elB.*omega_s;
Kmrr = Lm/Lrr;
R2 = Kmrr^2*Rr;
R1 = Rs+R2;
Ls_prime = Lss-Lm^2./Lrr;
Tr = Lrr/Rr;
k1 = omega_s.*Ls_prime/omega_el;
        v_dr = 0;                                           %  +.034 %for wr=0.95
        v_qr = 0;       %+0.00825; % for wr=1, -.035 %for wr=1.05, +.035 %for wr=0.95
        Vb = 1;
        v=10;
y0=[0.72367      0.33389      0.95043      0.11696       1.0042      1.0042 0 0]; 
[t01,x01] = ode15s(@open_loop_function, [0 20], y0);
y0=x01(end,:);
[t1,x1] = ode15s(@open_loop_function, [0:0.01:5], y0);
%disp(x1)
i_qs=x1(:,1);
i_ds=x1(:,2);
Ia=i_qs+1i*i_ds;
v_qs=real(Vb+Ia*j*Xe);
v_ds=imag(Vb+Ia*j*Xe);
Pe1=v_ds.*i_ds+v_qs.*i_qs;
Qe1=-v_qs.*i_ds+v_ds.*i_qs;

y0=x1(end,:);
v=10;
Vb = 0.4;
[t2,x2] = ode15s(@open_loop_function, [5:0.01:5.2], y0);
%disp(x2)
i_qs=x2(:,1);
i_ds=x2(:,2);
Ia=i_qs+1i*i_ds;
v_qs=real(Vb+Ia*j*Xe);
v_ds=imag(Vb+Ia*j*Xe);
Pe2=v_ds.*i_ds+v_qs.*i_qs;
Qe2=-v_qs.*i_ds+v_ds.*i_qs;

% plot([t1;t2],[x1(:,5);x2(:,5)],'b')
% plot([t1;t2],[Pe1;Pe2],'b')
% plot([t1;t2],[Qe1;Qe2],'b')

y0=x2(end,:);
v=10;
Vb = 0.85;

[t3,x3] = ode15s(@open_loop_function, [5.2:0.01:20], y0);

i_qs=x3(:,1);
i_ds=x3(:,2);
Ia=i_qs+1i*i_ds;
v_qs=real(Vb+Ia*j*Xe);
v_ds=imag(Vb+Ia*j*Xe);
Pe3=v_ds.*i_ds+v_qs.*i_qs;
Qe3=-v_qs.*i_ds+v_ds.*i_qs;

%hold on
plot([t1;t2;t3],[x1(:,5);x2(:,5);x3(:,5)],'b')
% plot([t1;t2;t3],[Qe1;Qe2;Qe3],'b')
% plot([t1;t2;t3],[Pe1;Pe2;Pe3],'b')

z1=[x1(:,5);x2(:,5);x3(:,5)];       %generator speed
z2=[Pe1;Pe2;Pe3];       % active power
z3=[Qe1;Qe2;Qe3];       %reactive power
z4=[x1(:,1);x2(:,1);x3(:,1)];
z5=[x1(:,2);x2(:,2);x3(:,2)];
z6=[x1(:,3);x2(:,3);x3(:,3)];
z7=[x1(:,4);x2(:,4);x3(:,4)];

% save z1 z1
% save z2 z2
% save z3 z3
% save z4 z4
% save z5 z5