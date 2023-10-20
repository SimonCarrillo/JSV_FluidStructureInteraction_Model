%% Function to get material constants prior to compatibility conditions per layer

function[du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z] = new_pagano_non(E,nu,i,j,hinit,hfinal,a,b,rho, omega)

n_thick = 101;

h = hfinal - hinit;

z = linspace(hinit,hfinal,n_thick)/h;

%%% Material Properties

lam = 1*nu/((1 + nu)*(1-2*nu));
mu = 1/(2*(1 + nu));

%%% Fourier Series Properties

p = i*pi*h/a;
q = j*pi*h/b;

a_bar = a/h;
tau = h*sqrt(rho/E);

omega_crit1 = sqrt(mu*(p^2+q^2));
omega_crit2 = sqrt((lam+2*mu)*(p^2+q^2));

%omega_crit1 = sqrt((mu*a_bar^4)*(p^2+q^2));
%omega_crit2 = sqrt((a_bar^4*(lam + 2*mu))*(p^2+q^2));

s1 = sqrt((p^2 + q^2- (omega^2)/(mu)));
s2 = sqrt((p^2 + q^2- (omega^2)/((lam + 2*mu))));

m1 = sqrt((omega^2/(mu) - p^2 - q^2));
m3 = sqrt((omega^2/((lam + 2*mu)) - p^2 - q^2));

%% non-dimensional for characteristic displacements and stresses

u_factor = 1;% 1/(a);
v_factor = 1;%1/(b);
w_factor = 1;%1/(hfinal);

xz_factor =1;%a^2/(100*E*hfinal^2);
yz_factor =1;%b^2/(100*E*hfinal^2);
zz_factor =1;%a^3/(100*E*hfinal^3);

%% case of k real

if (i == 0)
    
    if omega <= omega_crit1
        
        du1 = s1*exp(s1*z);
        du2 = s1*exp(-s1*z);
        du3 = 0*exp(s1*z);
        du4 = 0*exp(s1*z);
        du5 = 0*exp(s1*z);
        du6 = 0*exp(s1*z);
        
        du = [du1; du2; du3; du4; du5; du6]*u_factor;
        
        dv = zeros(6,n_thick);
        dw = zeros(6,n_thick);
        bzz = zeros(6,n_thick);
        byz = zeros(6,n_thick);
        bxx = zeros(6,n_thick);
        byy = zeros(6,n_thick);
        
        bxz1 = mu*s1^2*exp(s1*z);
        bxz2 = -mu*s1^2*exp(-s1*z);  
        bxz3 = 0*exp(s1*z);
        bxz4 = 0*exp(s1*z);
        bxz5 = 0*exp(s1*z);
        bxz6 = 0*exp(s1*z);
        
        bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6]*xz_factor;
        
        bxy1 = mu*q*s1*exp(s1*z);
        bxy2 = mu*q*s1*exp(-s1*z);  
        bxy3 = 0*exp(s1*z);
        bxy4 = 0*exp(s1*z);
        bxy5 = 0*exp(s1*z);
        bxy6 = 0*exp(s1*z);
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
        
    elseif omega > omega_crit1
        
        du1 = m1*cos(m1*z);
        du2 = m1*sin(m1*z);
        du3 = 0*cos(m1*z);
        du4 = 0*cos(m1*z);
        du5 = 0*cos(m1*z);
        du6 = 0*cos(m1*z);
        
        du = [du1; du2; du3; du4; du5; du6]*u_factor;
        
        dv = zeros(6,n_thick);
        dw = zeros(6,n_thick);
        bzz = zeros(6,n_thick);
        byz = zeros(6,n_thick);
         bxx = zeros(6,n_thick);
        byy = zeros(6,n_thick);
        
        bxz1 = -mu*m1^2*sin(m1*z);
       bxz2 = mu*m1^2*cos(m1*z);
        bxz3 = 0*cos(m1*z);
        bxz4 = 0*cos(m1*z);
        bxz5 = 0*cos(m1*z);
        bxz6 = 0*cos(m1*z);
        
        bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6]*xz_factor;
        
        
        bxy1 = mu*q*m1*cos(m1*z);
        bxy2 = mu*q*m1*sin(m1*z);  
        bxy3 = 0*exp(s1*z);
        bxy4 = 0*exp(s1*z);
        bxy5 = 0*exp(s1*z);
        bxy6 = 0*exp(s1*z);
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
    end
    
elseif (j == 0)
    
    if omega <= omega_crit1
        
        dv1 = s1*exp(s1*z);
        dv2 = s1*exp(-s1*z);
        dv3 = 0*exp(s1*z);
        dv4 = 0*exp(s1*z);
        dv5 = 0*exp(s1*z);
        dv6 = 0*exp(s1*z);
        
        dv = [dv1; dv2; dv3; dv4; dv5; dv6]*v_factor;
        
        du = zeros(6,n_thick);
        dw = zeros(6,n_thick);
        bzz = zeros(6,n_thick);
        bxz = zeros(6,n_thick);
         bxx = zeros(6,n_thick);
        byy = zeros(6,n_thick);
        
        byz1 = mu*s1^2*exp(s1*z);
        byz2 = -mu*s1^2*exp(-s1*z);
        byz3 = 0*exp(s1*z);
        byz4 = 0*exp(s1*z);
        byz5 = 0*exp(s1*z);
        byz6 = 0*exp(s1*z);
        
        byz= [byz1; byz2; byz3; byz4; byz5; byz6]*yz_factor;
        
        bxy1 = mu*p*s1*exp(s1*z);
        bxy2 = mu*p*s1*exp(-s1*z);  
        bxy3 = 0*exp(s1*z);
        bxy4 = 0*exp(s1*z);
        bxy5 = 0*exp(s1*z);
        bxy6 = 0*exp(s1*z);
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
    elseif  omega > omega_crit1
        
        dv1 = m1*cos(m1*z);
        dv2 = m1*sin(m1*z);
        dv3 = 0*cos(m1*z);
        dv4 = 0*cos(m1*z);
        dv5 = 0*cos(m1*z);
        dv6 = 0*cos(m1*z);
        
        dv = [dv1; dv2; dv3; dv4; dv5; dv6]*v_factor;
        
        du = zeros(6,n_thick);
        dw = zeros(6,n_thick);
        bzz = zeros(6,n_thick);
        bxz = zeros(6,n_thick);
         bxx = zeros(6,n_thick);
        byy = zeros(6,n_thick);
        
        byz1 = -mu*m1^2*sin(m1*z);
        byz2 = mu*m1^2*cos(m1*z);
        byz3 = 0*cos(m1*z);
        byz4 = 0*cos(m1*z);
        byz5 = 0*cos(m1*z);
        byz6 = 0*cos(m1*z);
        
        byz= [byz1; byz2; byz3; byz4; byz5; byz6]*yz_factor;
        
        bxy1 = mu*p*m1*cos(m1*z);
        bxy2 = mu*p*m1*sin(m1*z);  
        bxy3 = 0*exp(s1*z);
        bxy4 = 0*exp(s1*z);
        bxy5 = 0*exp(s1*z);
        bxy6 = 0*exp(s1*z);
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
    end
    
else
    
        if omega == 0
    
            c11 = 1*(1-nu)/((1-2*nu)*(1+nu));
    
            %%% constant definitions for displacement u
    
            du1 = exp(s1*z);
            du2 =  exp(-s1*z);
            du3 = z.*exp(s1*z);
            du4 =  z.*exp(-s1*z);
            du5 = 0*exp(s1*z);
            du6 = 0*exp(-s1*z);
    
            du = [du1; du2; du3; du4; du5; du6]*u_factor;
    
            %%% constant definitions for displacement v
    
            dv1 =  0*exp(s1*z);
            dv2 = 0*exp(s1*z);
            dv3 = q/p*z.*exp(s1*z);
            dv4 = q/p*z.*exp(-s1*z);
            dv5 = exp(s1*z);
            dv6 = exp(-s1*z);
    
            dv = [dv1; dv2; dv3; dv4; dv5; dv6]*v_factor;
    
            %%% constant definitions for displacement w
            dw1 = p*exp(s1*z)/s1;
            dw2 = -p*exp(-s1*z)/s1;
            dw3 = ((4*nu-3)/(p)+s1*z/p).*exp(s1*z);
            dw4 = ((4*nu-3)/(p)-s1*z/p).*exp(-s1*z);
            dw5 = q*exp(s1*z)/s1;
            dw6 = -q*exp(-s1*z)/s1;
    
    
            dw = [dw1; dw2; dw3; dw4; dw5; dw6]*w_factor;
    
    
            %%% constant definitions for traction zz
    
            bzz1 = c11*p*(1-2*nu)/(1-nu)*exp(s1*z);
            bzz2 = c11*p*(1-2*nu)/(1-nu)*exp(-s1*z);
            bzz3 =  c11*(s1*(1-2*nu))/(p*(1-nu))*exp(s1*z).*(s1*z-2*(1-nu));
            bzz4 = c11*(s1*(1-2*nu))/(p*(1-nu))*exp(-s1*z).*(s1*z+2*(1-nu));
            bzz5 =  c11*q*(1-2*nu)/(1-nu)*exp(s1*z);
            bzz6 =  c11*q*(1-2*nu)/(1-nu)*exp(-s1*z);
    
    
            bzz = [bzz1; bzz2; bzz3; bzz4; bzz5; bzz6]*zz_factor;
    
            %%% constant definitions for traction yz
    
            byz1 = c11*(1-2*nu)/(2*(1-nu))*p*q/s1*exp(s1*z);
            byz2 =  - c11*(1-2*nu)/(2*(1-nu))*p*q/s1*exp(-s1*z);
            byz3 =   c11*(1-2*nu)/(1-nu)*(2*nu-1+s1*z).*q/p.*exp(s1*z);
            byz4 =    c11*(1-2*nu)/(1-nu)*(2*nu-1-s1*z).*q/p.*exp(-s1*z);
            byz5 =    c11*(1-2*nu)/(2*(1-nu))*(q^2/s1 + s1)*exp(s1*z);
            byz6 =   - c11*(1-2*nu)/(2*(1-nu))*(q^2/s1 + s1)*exp(-s1*z);
    
    
            byz = [byz1; byz2; byz3; byz4; byz5; byz6]*yz_factor;
    
            %%% constant definitions for traction xz
    
            bxz1 = c11*(1-2*nu)/(2*(1-nu))*(p^2/s1 + s1)*exp(s1*z);
            bxz2 = - c11*(1-2*nu)/(2*(1-nu))*(p^2/s1 + s1)*exp(-s1*z);
            bxz3 =  c11*(1-2*nu)/(1-nu)*(2*nu-1+s1*z).*exp(s1*z);
            bxz4 =    c11*(1-2*nu)/(1-nu)*(2*nu-1-s1*z).*exp(-s1*z);
            bxz5 =  c11*(1-2*nu)/(2*(1-nu))*p*q/s1*exp(s1*z);
            bxz6 =  - c11*(1-2*nu)/(2*(1-nu))*p*q/s1*exp(-s1*z);
    
    
            bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6]*xz_factor;
    
    
    elseif (omega <= omega_crit1)
        
        %%% constant definitions for displacement u
        
        du1 = s1*exp(s1*z);
        du2 =  0*exp(s1*z);
        du3 = s1*exp(-s1*z);
        du4 =  0*exp(-s1*z);
        du5 = p*exp(s2*z);
        du6 = p*exp(-s2*z);
        
        du = [du1; du2; du3; du4; du5; du6];
        
        %%% constant definitions for displacement v
        
        dv1 =  0*exp(s1*z);
        dv2 = s1*exp(s1*z);
        dv3 = 0*exp(-s1*z);
        dv4 = s1*exp(-s1*z);
        dv5 = q*exp(s2*z);
        dv6 = q*exp(-s2*z);
        
        dv = [dv1; dv2; dv3; dv4; dv5; dv6];
        
        %%% constant definitions for displacement w
        dw1 = p*exp(s1*z);
        dw2 = q*exp(s1*z);
        dw3 = -p*exp(-s1*z);
        dw4 = -q*exp(-s1*z);
        dw5 = s2*exp(s2*z);
        dw6 = -s2*exp(-s2*z);
        
        
        dw = [dw1; dw2; dw3; dw4; dw5; dw6];
        
        
        %%% constant definitions for traction zz
        
        bzz1 = 2*mu*p*s1*exp(s1*z);
        bzz2 = 2*mu*q*s1*exp(s1*z);
        bzz3 =  2*mu*p*s1*exp(-s1*z);
        bzz4 =  2*mu*q*s1*exp(-s1*z);
        bzz5 =  (-lam*(p^2 + q^2) +(lam + 2*mu)*s2^2)*exp(s2*z);
        bzz6 =  (-lam*(p^2 + q^2) +(lam + 2*mu)*s2^2)*exp(-s2*z);
        
        
        bzz = [bzz1; bzz2; bzz3; bzz4; bzz5; bzz6];
        
        %%% constant definitions for traction yz
        
        byz1 = (mu*p*q)*exp(s1*z);
        byz2 =    (mu*(q^2+s1^2))*exp(s1*z);
        byz3 =    (-mu*p*q)*exp(-s1*z);
        byz4 =    (-mu*(q^2+s1^2))*exp(-s1*z);
        byz5 =    (2*mu*q*s2)*exp(s2*z);
        byz6 =    (-2*mu*q*s2)*exp(-s2*z);
        
        
        byz = [byz1; byz2; byz3; byz4; byz5; byz6];
        
        %%% constant definitions for traction xz
        
        bxz1 = mu*(p^2+s1^2)*exp(s1*z);
        bxz2 =    mu*p*q*exp(s1*z);
        bxz3 =    -mu*(p^2+s1^2)*exp(-s1*z);
        bxz4 =    -mu*p*q*exp(-s1*z);
        bxz5 =   2*mu*p*s2*exp(s2*z);
        bxz6 =    -2*mu*p*s2*exp(-s2*z);
        
        
        bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6];
        
                %%% constant definitions for traction xx
        
        bxx1 = -2*mu*p*s1*exp(s1*z);
        bxx2 =    0*exp(s1*z);
        bxx3 =    -2*mu*p*s1*exp(-s1*z);
        bxx4 =    0*exp(-s1*z);
        bxx5 =   (-(lam+2*mu)*p^2-lam*q^2+lam*s2^2)*exp(s2*z);
        bxx6 =    (-(lam+2*mu)*p^2-lam*q^2+lam*s2^2)*exp(-s2*z);
      
        
        bxx= [bxx1; bxx2; bxx3; bxx4; bxx5; bxx6];
        
                        %%% constant definitions for traction yy
        
        byy1 = 0*exp(s1*z);
        byy2 =    -2*mu*q*s1*exp(s1*z);
        byy3 =    0*exp(-s1*z);
        byy4 =    -2*mu*q*s1*exp(-s1*z);
        byy5 =   (-(lam+2*mu)*q^2-lam*p^2+lam*s2^2)*exp(s2*z);
        byy6 =    (-(lam+2*mu)*q^2-lam*p^2+lam*s2^2)*exp(-s2*z);
      
        
        byy= [byy1; byy2; byy3; byy4; byy5; byy6];
        
                                %%% constant definitions for traction xy
        
        bxy1 = mu*q*s1*exp(s1*z);
        bxy2 =    mu*p*s1*exp(s1*z);
        bxy3 =    mu*q*s1*exp(-s1*z);
        bxy4 =    mu*p*s1*exp(-s1*z);
        bxy5 =   2*mu*p*q*exp(s2*z);
        bxy6 =   2*mu*p*q*exp(-s2*z);
      
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
        
        %% case of k1 imag and k3 real
        
    elseif (omega > omega_crit1) && (omega <= omega_crit2)
        
        
        du1 = m1*cos(m1*z);
        du2 =  0*cos(m1*z);
        du3 = m1*sin(m1*z);
        du4 =  0*sin(m1*z);
        du5 = p*exp(s2*z);
        du6 = p*exp(-s2*z);
        
        du = [du1; du2; du3; du4; du5; du6];
        
        %%% constant definitions for displacement v
        
        dv1 =  0*cos(m1*z);
        dv2 = m1*cos(m1*z);
        dv3 = 0*sin(m1*z);
        dv4 = m1*sin(m1*z);
        dv5 = q*exp(s2*z);
        dv6 = q*exp(-s2*z);
        
        dv = [dv1; dv2; dv3; dv4; dv5; dv6];
        
        %%% constant definitions for displacement w
        dw1 = p*sin(m1*z);
        dw2 = q*sin(m1*z);
        dw3 = -p*cos(m1*z);
        dw4 = -q*cos(m1*z);
        dw5 = s2*exp(s2*z);
        dw6 = -s2*exp(-s2*z);
        
        
        dw = [dw1; dw2; dw3; dw4; dw5; dw6];
        
        
        %%% constant definitions for traction zz
        
        bzz1 = 2*mu*p*m1*cos(m1*z);
        bzz2 = 2*mu*q*m1*cos(m1*z);
        bzz3 =  2*mu*p*m1*sin(m1*z);
        bzz4 =  2*mu*q*m1*sin(m1*z);
        bzz5 =  (-lam*(p^2 + q^2) +(lam + 2*mu)*s2^2)*exp(s2*z);
        bzz6 =  (-lam*(p^2 + q^2) +(lam + 2*mu)*s2^2)*exp(-s2*z);
        
        
        bzz = [bzz1; bzz2; bzz3; bzz4; bzz5; bzz6];
        
        %%% constant definitions for traction yz
        
        byz1 = (mu*p*q)*sin(m1*z);
        byz2 =    (mu*(q^2-m1^2))*sin(m1*z);
        byz3 =    (-mu*p*q)*cos(m1*z);
        byz4 =    (mu*(m1^2-q^2))*cos(m1*z);
        byz5 =    (2*mu*q*s2)*exp(s2*z);
        byz6 =    (-2*mu*q*s2)*exp(-s2*z);
        
        
        byz = [byz1; byz2; byz3; byz4; byz5; byz6];
        
        %%% constant definitions for traction xz
        
        bxz1 = mu*(p^2-m1^2)*sin(m1*z);
        bxz2 =    mu*p*q*sin(m1*z);
        bxz3 =    mu*(m1^2-p^2)*cos(m1*z);
        bxz4 =    -mu*p*q*cos(m1*z);
        bxz5 =   2*mu*p*s2*exp(s2*z);
        bxz6 =    -2*mu*p*s2*exp(-s2*z);
        
        
        bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6];
        
                               %%% constant definitions for traction xx
        
        bxx1 = -2*mu*p*m1*cos(m1*z);
        bxx2 =    0*cos(m1*z);
        bxx3 =    -2*mu*p*m1*sin(m1*z);
        bxx4 =    0*sin(m1*z);
        bxx5 =   (-(lam+2*mu)*p^2-lam*q^2+lam*s2^2)*exp(s2*z);
        bxx6 =    (-(lam+2*mu)*p^2-lam*q^2+lam*s2^2)*exp(-s2*z);
      
        
        bxx= [bxx1; bxx2; bxx3; bxx4; bxx5; bxx6];
        
                        %%% constant definitions for traction yy
        
        byy1 = 0*cos(m1*z);
        byy2 =    -2*mu*q*m1*cos(m1*z);
        byy3 =    0*sin(m1*z);
        byy4 =    -2*mu*q*m1*sin(m1*z);
        byy5 =   (-(lam+2*mu)*q^2-lam*p^2+lam*s2^2)*exp(s2*z);
        byy6 =    (-(lam+2*mu)*q^2-lam*p^2+lam*s2^2)*exp(-s2*z);
      
        
        byy= [byy1; byy2; byy3; byy4; byy5; byy6];
        
                                %%% constant definitions for traction xy
        
        bxy1 = mu*q*m1*cos(m1*z);
        bxy2 =    mu*p*m1*cos(m1*z);
        bxy3 =    mu*q*m1*sin(m1*z);
        bxy4 =    mu*p*m1*sin(m1*z);
         bxy5 =   2*mu*p*q*exp(s2*z);
        bxy6 =   2*mu*p*q*exp(-s2*z);
      
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
        
        
        
        %% case of k imaginary
        
    elseif omega > omega_crit2
        
        du1 = m1*cos(m1*z);
        du2 =  0*cos(m1*z);
        du3 = m1*sin(m1*z);
        du4 =  0*sin(m1*z);
        du5 = p*cos(m3*z);
        du6 = p*sin(m3*z);
        
        du = [du1; du2; du3; du4; du5; du6];
        
        %%% constant definitions for displacement v
        
        dv1 =  0*cos(m1*z);
        dv2 = m1*cos(m1*z);
        dv3 = 0*sin(m1*z);
        dv4 = m1*sin(m1*z);
        dv5 = q*cos(m3*z);
        dv6 = q*sin(m3*z);
        
        dv = [dv1; dv2; dv3; dv4; dv5; dv6];
        
        %%% constant definitions for displacement w
        dw1 = p*sin(m1*z);
        dw2 = q*sin(m1*z);
        dw3 = -p*cos(m1*z);
        dw4 = -q*cos(m1*z);
        dw5 = -m3*sin(m3*z);
        dw6 = m3*cos(m3*z);
        
        
        dw = [dw1; dw2; dw3; dw4; dw5; dw6];
        
        
        %%% constant definitions for traction zz
        
        bzz1 = 2*mu*p*m1*cos(m1*z);
        bzz2 = 2*mu*q*m1*cos(m1*z);
        bzz3 =  2*mu*p*m1*sin(m1*z);
        bzz4 =  2*mu*q*m1*sin(m1*z);
        bzz5 =  -(lam*(p^2 + q^2) +(lam + 2*mu)*m3^2)*cos(m3*z);
        bzz6 =  -(lam*(p^2 + q^2) +(lam + 2*mu)*m3^2)*sin(m3*z);
        
        
        bzz = [bzz1; bzz2; bzz3; bzz4; bzz5; bzz6];
        
        %%% constant definitions for traction yz
        
        byz1 = (mu*p*q)*sin(m1*z);
        byz2 =    (mu*(q^2-m1^2))*sin(m1*z);
        byz3 =    (-mu*p*q)*cos(m1*z);
        byz4 =    (mu*(m1^2-q^2))*cos(m1*z);
        byz5 =    -(2*mu*q*m3)*sin(m3*z);
        byz6 =    (2*mu*q*m3)*cos(m3*z);
        
        
        byz = [byz1; byz2; byz3; byz4; byz5; byz6];
        
        %%% constant definitions for traction xz
        
        bxz1 = mu*(p^2-m1^2)*sin(m1*z);
        bxz2 =    mu*p*q*sin(m1*z);
        bxz3 =    mu*(m1^2-p^2)*cos(m1*z);
        bxz4 =    -mu*p*q*cos(m1*z);
        bxz5 =   -2*mu*p*m3*sin(m3*z);
        bxz6 =    2*mu*p*m3*cos(m3*z);
        
        
        bxz= [bxz1; bxz2; bxz3; bxz4; bxz5; bxz6];
        
        
                       %%% constant definitions for traction xx
        
        bxx1 = -2*mu*p*m1*cos(m1*z);
        bxx2 =    0*cos(m1*z);
        bxx3 =    -2*mu*p*m1*sin(m1*z);
        bxx4 =    0*sin(m1*z);
        bxx5 =   (-(lam+2*mu)*p^2-lam*q^2+lam*m3^2)*cos(m3*z);
        bxx6 =    (-(lam+2*mu)*p^2-lam*q^2+lam*m3^2)*sin(m3*z);
      
        
        bxx= [bxx1; bxx2; bxx3; bxx4; bxx5; bxx6];
        
                        %%% constant definitions for traction yy
        
        byy1 = 0*cos(m1*z);
        byy2 =    -2*mu*q*m1*cos(m1*z);
        byy3 =    0*sin(m1*z);
        byy4 =    -2*mu*q*m1*sin(m1*z);
        byy5 =   (-(lam+2*mu)*q^2-lam*p^2+lam*m3^2)*cos(m3*z);
        byy6 =    (-(lam+2*mu)*q^2-lam*p^2+lam*m3^2)*sin(m3*z);
      
        
        byy= [byy1; byy2; byy3; byy4; byy5; byy6];
        
                                %%% constant definitions for traction xy
        
        bxy1 = mu*q*m1*cos(m1*z);
        bxy2 =    mu*p*m1*cos(m1*z);
        bxy3 =    mu*q*m1*sin(m1*z);
        bxy4 =    mu*p*m1*sin(m1*z);
        bxy5 =   2*mu*p*q*cos(m3*z);
        bxy6 =   2*mu*p*q*sin(m3*z);
      
        
        bxy= [bxy1; bxy2; bxy3; bxy4; bxy5; bxy6];
        
        
    end
    
end

z= z'; % neccesary to redefine the total z coordinate outside

end