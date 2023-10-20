clc
clear all

%%% 3D Fluid-Structure Interaction Aluminum Panel
%%% NYU Tandon School of Engineering
%%% Author: Simon Carrillo Segura

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FSI = 1; % 0 for dry and 1 for wet

flag_save_dry = 0;  % 1 for saving results
flag_save_wet = 1;  % 1 for saving results

a_test = 1;
H_test = 2;

number_i= 10; % number i,j
number_j= 10; % number i,j

number_m= 10; % number of m,n
number_n= 10; % number of m,n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Do not touch below this line %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_parametric = [0.04 0.1 0.2 0.4];
H_parametric = [0.04 0.1 0.2 0.4];

rho_fluid = 1000; %1000

% To change for the parametric analysis
H = H_parametric(H_test); % tank height

a =  a_parametric(a_test); % in meters plate dimensions

b = a;
hf = 0.02;

Ef= 69*10^9; % in Pa
nuf= 0.3;
rhof = 2.702*10^3; % in kg/m3

omega_factor = sqrt(rhof/Ef)*hf;

n_z= 6; % number of unknowns in the z-direction (6 unknowns per layer)
n_thick = 101; % number of datapoints in the z-direction for plotting %%% 100000*2*pi*omega_factor

omega = linspace(0.001, 90000*2*pi*omega_factor, 3000); %90k, 40k, 17k, 6.5k %omega range to look into

%% Vaccuum solution
if FSI == 0
    
    [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_n, dummy_eig,  dummy_dia,condition_A] = Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega, number_i,number_j, number_m, number_n, n_thick, n_z);
    
else
    %% Fluid-Vaccum solution
    
    [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_n, dummy_eig,  dummy_dia,condition_A] = Fluid_Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega, rho_fluid, H, number_i,number_j, number_m, number_n, n_thick, n_z);
    
end
%% Figures d_n vs omega

figure()
semilogy(omega/(2*pi)/omega_factor, condition_A,'LineWidth',2)
xlabel('$f$ $[Hz]$','interpreter','latex','FontSize',14,'LineWidth',2)
ylabel('condition of A','interpreter','latex','FontSize',14)
hold on
grid on

%% Look for local minima and input to Secant Method and Local Minima from d_n values

dn_abs = abs(d_n);

TF_n = islocalmin(dn_abs);
TF_con = islocalmax(condition_A);

% Obtain all local minima for possible eigenvalues

for i=1:length(omega)
    
    if (TF_n(i) ==1)
        
        omega_min(i) = omega(i);
        
    end
    
end

omegas = nonzeros(omega_min);
omegas_Hz = nonzeros(omega_min)/(2*pi)/omega_factor;

%% Refine approximation of eigenfrequencies

for iter_omegas=1:length(omegas)
    
    omega0 = omegas(iter_omegas)-0.002; % 0.0001 for thin plates
    omega1 = omegas(iter_omegas)+0.002; % 0.0001 for thin plates
    
    omegas_refined = linspace(omega0,omega1, 100);
    
    if FSI == 0
        
        [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_n1, dummy_eig1,  dummy_dia1,condition_A1] = Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omegas_refined, number_i,number_j, number_m, number_n, n_thick, n_z);
        
    else
        
        [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_n1, dummy_eig1,  dummy_dia1,condition_A1] = Fluid_Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omegas_refined, rho_fluid, H, number_i,number_j, number_m, number_n, n_thick, n_z);
        
    end
    
    TF_1 = islocalmin(abs(d_n1));
    
    %% Obtain all local minima for possible eigenvalues
    
    for iter=1:length(omegas_refined)
        
        if (TF_1(iter) ==1) 
            
            omega_final(iter) = omegas_refined(iter);
            dn_final(iter) = d_n1(iter);
            
        else
            
            omega_final(iter) =0;
            dn_final(iter) = 0;
            
        end
        
    end
    
    omega_final = nonzeros(omega_final)
    dn_final = nonzeros(dn_final);
    
    figure()
    plot(omegas_refined, d_n1)

    
    %% Mode shape reconstruction
    
    if FSI == 0
        
        [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_nfinal, dummy_eigfinal(:,:,iter_omegas),  dummy_diafinal(:,iter_omegas),condition_Afinal] = Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega_final(end), number_i,number_j, number_m, number_n, n_thick, n_z);
        
    else
        
        [du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, A, d_nfinal, dummy_eigfinal(:,:,iter_omegas),  dummy_diafinal(:,iter_omegas),condition_Afinal] = Fluid_Vacumm_plate_nonzero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega_final(end), rho_fluid, H, number_i,number_j, number_m, number_n, n_thick, n_z);

    end
    
    
    %% analyze if we have a multiple root
    
    x= linspace(0,a,101)/hf;
    y = linspace(0,b,101)/hf;
    
    z_new = z(1,:);
    
    a_bar = a/hf;
    b_bar = b/hf;
    
    u_int = nan(number_i*number_j,length(z_new));
    v_int = nan(number_i*number_j,length(z_new));
    w_int = nan(number_i*number_j,length(z_new));
    sigmazz_int = nan(number_i*number_j,length(z_new));
    sigmaxz_int = nan(number_i*number_j,length(z_new));
    sigmayz_int = nan(number_i*number_j,length(z_new));
    sigmaxx_int = nan(number_i*number_j,length(z_new));
    sigmayy_int = nan(number_i*number_j,length(z_new));
    sigmaxy_int = nan(number_i*number_j,length(z_new));
    
    U_total = nan(length(x),length(y),length(z_new));
    V_total = nan(length(x),length(y),length(z_new));
    W_total = nan(length(x),length(y),length(z_new));
    sigmazz_total = nan(length(x),length(y),length(z_new));
    sigmaxz_total = nan(length(x),length(y),length(z_new));
    sigmayz_total = nan(length(x),length(y),length(z_new));
    sigmaxx_total = nan(length(x),length(y),length(z_new));
    sigmayy_total = nan(length(x),length(y),length(z_new));
    sigmaxy_total = nan(length(x),length(y),length(z_new));
    
    u_int2 = nan(number_i*number_j,length(z_new));
    v_int2 = nan(number_i*number_j,length(z_new));
    w_int2 = nan(number_i*number_j,length(z_new));
    sigmazz_int2 = nan(number_i*number_j,length(z_new));
    sigmaxz_int2 = nan(number_i*number_j,length(z_new));
    sigmayz_int2 = nan(number_i*number_j,length(z_new));
    sigmaxx_int2 = nan(number_i*number_j,length(z_new));
    sigmayy_int2 = nan(number_i*number_j,length(z_new));
    sigmaxy_int2 = nan(number_i*number_j,length(z_new));
    
    U_total2 = nan(length(x),length(y),length(z_new));
    V_total2 = nan(length(x),length(y),length(z_new));
    W_total2 = nan(length(x),length(y),length(z_new));
    sigmazz_total2 = nan(length(x),length(y),length(z_new));
    sigmaxz_total2 = nan(length(x),length(y),length(z_new));
    sigmayz_total2 = nan(length(x),length(y),length(z_new));
    sigmaxx_total2 = nan(length(x),length(y),length(z_new));
    sigmayy_total2 = nan(length(x),length(y),length(z_new));
    sigmaxy_total2 = nan(length(x),length(y),length(z_new));
    
    eps = 10^-9;
    
    if abs(dummy_diafinal(end,iter_omegas) - dummy_diafinal(end-1,iter_omegas))< eps % case of multiple roots
        
        U_root = dummy_eigfinal(:,end,iter_omegas);
        U_root = U_root/max(abs(U_root));
        
        U_root2 = dummy_eigfinal(:,end-1, iter_omegas);
        U_root2 = U_root2/max(abs(U_root2));
        
        display('multiple root')
        
        %% computation of displacement and stresses in the z-direction
        
        for l= 1:length(z_new)
            
            for r= 1: (number_i*number_j)
                
                u_int(r,l) = du(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                v_int(r,l) = dv(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                w_int(r,l) = dw(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmazz_int(r,l) = bzz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxz_int(r,l) = bxz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmayz_int(r,l) = byz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxx_int(r,l) = bxx(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmayy_int(r,l) = byy(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxy_int(r,l) = bxy(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                
            end
            
        end
        
        %% computation of global displacement and stresses from the contribution of all i,j
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    U_total(t,o,l) = 0;
                    V_total(t,o,l) = 0;
                    W_total(t,o,l) = 0;
                    sigmazz_total(t,o,l) =0;
                    sigmaxz_total(t,o,l) =0;
                    sigmayz_total(t,o,l) =0;
                    sigmaxx_total(t,o,l) =0;
                    sigmayy_total(t,o,l) =0;
                    sigmaxy_total(t,o,l) =0;
                    
                    for i = 1:number_i
                        
                        for j=1:number_j
                            
                            r = (i-1)*number_i+j;
                            
                            U_total(t,o,l) = U_total(t,o,l) + u_int(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            V_total(t,o,l) = V_total(t,o,l) + v_int(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            W_total(t,o,l) = W_total(t,o,l) + w_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmazz_total(t,o,l) = sigmazz_total(t,o,l) + sigmazz_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxz_total(t,o,l) = sigmaxz_total(t,o,l) + sigmaxz_int(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayz_total(t,o,l) = sigmayz_total(t,o,l) + sigmayz_int(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            sigmaxx_total(t,o,l) = sigmaxx_total(t,o,l) + sigmaxx_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayy_total(t,o,l) = sigmayy_total(t,o,l) + sigmayy_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxy_total(t,o,l) = sigmaxy_total(t,o,l) + sigmaxy_int(r,l).*cos((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);

                        end
                    end
                    
                end
                
            end
            
        end
        
        %% computation of displacement and stresses in the z-direction
        
        
        for l= 1:length(z_new)
            
            for r= 1: (number_i*number_j)
                
                u_int2(r,l) = du(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                v_int2(r,l) = dv(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                w_int2(r,l) = dw(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmazz_int2(r,l) = bzz(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmaxz_int2(r,l) = bxz(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmayz_int2(r,l) = byz(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmaxx_int2(r,l) = bxx(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmayy_int2(r,l) = byy(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                sigmaxy_int2(r,l) = bxy(r,:,l)*U_root2(((r-1)*n_z+1):(r*n_z));
                
            end
            
        end
        
        %% computation of global displacement and stresses from the contribution of all i,j
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    U_total2(t,o,l) = 0;
                    V_total2(t,o,l) = 0;
                    W_total2(t,o,l) = 0;
                    sigmazz_total2(t,o,l) =0;
                    sigmaxz_total2(t,o,l) =0;
                    sigmayz_total2(t,o,l) =0;
                    sigmaxx_total2(t,o,l) =0;
                    sigmayy_total2(t,o,l) =0;
                    sigmaxy_total2(t,o,l) =0;
                    
                    for i = 1:number_i
                        
                        for j=1:number_j
                            
                            r = (i-1)*number_i+j;
                            
                            U_total2(t,o,l) = U_total2(t,o,l) + u_int2(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            V_total2(t,o,l) = V_total2(t,o,l) + v_int2(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            W_total2(t,o,l) = W_total2(t,o,l) + w_int2(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmazz_total2(t,o,l) = sigmazz_total2(t,o,l) + sigmazz_int2(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxz_total2(t,o,l) = sigmaxz_total2(t,o,l) + sigmaxz_int2(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayz_total2(t,o,l) = sigmayz_total2(t,o,l) + sigmayz_int2(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            sigmaxx_total2(t,o,l) = sigmaxx_total2(t,o,l) + sigmaxx_int2(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayy_total2(t,o,l) = sigmayy_total2(t,o,l) + sigmayy_int2(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxy_total2(t,o,l) = sigmaxy_total2(t,o,l) + sigmaxy_int2(r,l).*cos((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);

                        end
                    end
                    
                end
                
            end
            
        end
        
                        %% Through thickness figures
            
                w = reshape(W_total(50,50,:),[1,n_thick]);
                sigzz = reshape(sigmazz_total(50,50,:),[1,n_thick]);
            
                u = reshape(U_total(1,50,:),[1,n_thick]);
                sigxz = reshape(sigmaxz_total(1,50,:),[1,n_thick]);

        %%%%%%%%
        %% Save files
        if flag_save_dry > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            name2 = "D";
            test = strcat(name2,ch2,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'u', 'w', 'sigxz','sigzz',  'U_total','V_total','W_total', 'sigmazz_total','sigmaxz_total', 'sigmayz_total','sigmaxx_total','sigmayy_total','sigmaxy_total', 'W_total2', 'U_total2','V_total2', 'sigmazz_total2','sigmaxz_total2','sigmayz_total2', 'sigmaxx_total2','sigmayy_total2','sigmaxy_total2', 'du', 'dv', 'dw', 'bzz', 'bxz', 'byz', 'bxx','byy','bxy', 'U_root','U_root2');
            
        end
        
        if flag_save_wet > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            ch3 = int2str(H_test);
            %name2 = "wetmodes";
            name2 = "W";
            test = strcat(name2,ch2,ch3,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'u', 'w', 'sigxz','sigzz',  'U_total','V_total','W_total', 'sigmazz_total','sigmaxz_total', 'sigmayz_total','sigmaxx_total','sigmayy_total','sigmaxy_total', 'W_total2', 'U_total2','V_total2', 'sigmazz_total2','sigmaxz_total2','sigmayz_total2', 'sigmaxx_total2','sigmayy_total2','sigmaxy_total2', 'du', 'dv', 'dw', 'bzz', 'bxz', 'byz', 'bxx','byy','bxy', 'U_root','U_root2');
            
        end
        
        %%%%%%%%%%%%%%%%%%%
        
    else
        
          U_root = dummy_eigfinal(:,end,iter_omegas);
            
            U_root = U_root/max(abs(U_root));
            display('simple root')
            
        %% computation of displacement and stresses in the z-direction
        
        for l= 1:length(z_new)
            
            for r= 1: (number_i*number_j)
                
                u_int(r,l) = du(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                v_int(r,l) = dv(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                w_int(r,l) = dw(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmazz_int(r,l) = bzz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxz_int(r,l) = bxz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmayz_int(r,l) = byz(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxx_int(r,l) = bxx(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmayy_int(r,l) = byy(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                sigmaxy_int(r,l) = bxy(r,:,l)*U_root(((r-1)*n_z+1):(r*n_z));
                
            end
            
        end
        
        %% computation of global displacement and stresses from the contribution of all i,j
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    U_total(t,o,l) = 0;
                    V_total(t,o,l) = 0;
                    W_total(t,o,l) = 0;
                    sigmazz_total(t,o,l) =0;
                    sigmaxz_total(t,o,l) =0;
                    sigmayz_total(t,o,l) =0;
                    sigmaxx_total(t,o,l) =0;
                    sigmayy_total(t,o,l) =0;
                    sigmaxy_total(t,o,l) =0;
                    
                    for i = 1:number_i
                        
                        for j=1:number_j
                            
                            r = (i-1)*number_i+j;
                            
                            U_total(t,o,l) = U_total(t,o,l) + u_int(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            V_total(t,o,l) = V_total(t,o,l) + v_int(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            W_total(t,o,l) = W_total(t,o,l) + w_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmazz_total(t,o,l) = sigmazz_total(t,o,l) + sigmazz_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxz_total(t,o,l) = sigmaxz_total(t,o,l) + sigmaxz_int(r,l).*cos((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayz_total(t,o,l) = sigmayz_total(t,o,l) + sigmayz_int(r,l).*sin((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);
                            sigmaxx_total(t,o,l) = sigmaxx_total(t,o,l) + sigmaxx_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmayy_total(t,o,l) = sigmayy_total(t,o,l) + sigmayy_int(r,l).*sin((i)*pi*x(t)/a_bar).*sin((j)*pi*y(o)/b_bar);
                            sigmaxy_total(t,o,l) = sigmaxy_total(t,o,l) + sigmaxy_int(r,l).*cos((i)*pi*x(t)/a_bar).*cos((j)*pi*y(o)/b_bar);

                        end
                    end
                    
                end
                
            end
            
        end            
       %% Through thickness figures
            
                w = reshape(W_total(50,50,:),[1,n_thick]);
                sigzz = reshape(sigmazz_total(50,50,:),[1,n_thick]);
            
                u = reshape(U_total(1,50,:),[1,n_thick]);
                sigxz = reshape(sigmaxz_total(1,50,:),[1,n_thick]);
            
        %% Save files
        if flag_save_dry > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            name2 = "D";
            test = strcat(name2,ch2,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'u', 'w', 'sigxz','sigzz', 'W_total', 'U_total','V_total', 'sigmazz_total','sigmaxz_total','sigmayz_total','sigmaxx_total','sigmayy_total','sigmaxy_total', 'du', 'dv', 'dw', 'bzz', 'bxz', 'byz','bxx','byy','bxy', 'U_root');
            
        end
        
        if flag_save_wet > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            ch3 = int2str(H_test);
            %name2 = "wetmodes";
            name2 = "W";
            test = strcat(name2,ch2,ch3,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'u', 'w', 'sigxz','sigzz', 'W_total', 'U_total','V_total', 'sigmazz_total','sigmaxz_total','sigmayz_total','sigmaxx_total','sigmayy_total','sigmaxy_total', 'du', 'dv', 'dw', 'bzz', 'bxz', 'byz','bxx','byy','bxy', 'U_root');
            
        end
        
    end
    
    clear omega_final
    
end        
