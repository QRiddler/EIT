
bk=0.1;
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

