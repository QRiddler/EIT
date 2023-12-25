function dx = open_loop_function(t,x)
global Xe e_prime_qs e_prime_ds i_qs i_ds omega_r thetatw omega_t Vs gamae Qtogrid v_qs v_ds v_dr  v_qr  
global v Vb omega_s omega_elB omega_el Kmrr R2 R1 Ls_prime Tr Pt k c Ht Hg Xm Rs Rr Lm Lss Lrr c bk

i_qs = x(1);
i_ds = x(2);
e_prime_qs = x(3);
e_prime_ds = x(4);
omega_r = x(5);
omega_t = x(6);
theta_r=x(7);
theta_t=x(8);

dx = zeros(8,1);
k1=omega_s*Ls_prime/omega_el;
Ia=i_qs+1i*i_ds;
dx(1)=(-R1*i_qs+omega_s.*Ls_prime.*i_ds+(omega_r./omega_s).*e_prime_qs-(1/(Tr.*omega_s)).*e_prime_ds-real(Vb+Ia*j*Xe)+Kmrr.*v_qr)/k1;
dx(2)=(-omega_s.*Ls_prime.*i_qs-R1*i_ds+(1/(Tr.*omega_s)).*e_prime_qs+(omega_r./omega_s).*e_prime_ds-(imag(Vb+Ia*j*Xe))+Kmrr.*v_dr)/k1;
dx(3)=(R2*i_ds-(1/(Tr.*omega_s)).*e_prime_qs+(1-(omega_r./omega_s)).*e_prime_ds-Kmrr.*v_dr).*omega_el;
dx(4)=(-R2*i_qs-(1-(omega_r./omega_s)).*e_prime_qs-(1/(Tr.*omega_s)).*e_prime_ds+Kmrr.*v_qr).*omega_el;
        
        Te = (e_prime_qs.*i_qs./omega_s)+(e_prime_ds.*i_ds./omega_s);
        
% Backlash model
theta_d=theta_t-theta_r;
if (theta_d+(c/k)*(omega_t-omega_r)*omega_elB)>bk
    Tsh = k*(theta_d-bk)+c.*(omega_t-omega_r).*omega_elB;
elseif (theta_d+(c/k)*(omega_t-omega_r)*omega_elB)<-bk
    Tsh=k*(theta_d+bk)+c.*(omega_t-omega_r).*omega_elB;
elseif abs(theta_d+(c/k)*(omega_t-omega_r)*omega_elB)<=bk
    Tsh=0;
end
%% Dead Zone
% 
% theta_d=theta_t-theta_r;
% if theta_d>bk
%     Tsh = k*(theta_d-bk)+c.*(omega_t-omega_r).*omega_elB;
% elseif theta_d<-bk
%     Tsh=k*(theta_d+bk)+c.*(omega_t-omega_r).*omega_elB;
% elseif abs(theta_d)<=bk
%     Tsh=0;
% end
        
dx(5) = (Tsh-Te)./(2*Hg);

beta = 0;
kopt = 1;
wr = omega_r*4.5;
lembda = (wr*(43.3/2)./v);
lembdai = 1./(1./(lembda+.08*beta)-.035/(beta^3+1));
Cp = (0.5176*(116./lembdai-0.4*beta-5).*exp(-21./lembdai)+0.0068.*lembda);
cpmax = 0.48;
cppu = (1/cpmax)*Cp;
vpu = v/12;
Pt = kopt*cppu*(vpu)^3;

dx(6) = ((Pt./omega_t)-Tsh)./(2*Ht);
dx(7) = (omega_r-omega_s).*omega_elB;
dx(8) = (omega_t-omega_s).*omega_elB;
end

