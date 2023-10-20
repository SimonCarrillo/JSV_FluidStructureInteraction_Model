clc
clear all

%%% 3D Fluid-Structure Interaction Aluminum Panel
%%% NYU Tandon School of Engineering
%%% Author: Simon Carrillo Segura

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_save_dry = 1;  % 1 for saving results
flag_save_wet = 0;  % 1 for saving results

a_test = 1;
H_test = 1;

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
n_thick = 101; % number of datapoints in the z-direction for plotting

omega = linspace(10*2*pi*omega_factor, 90000*2*pi*omega_factor, 1000); %110k, 40k, 17k, 6.5k %omega range to look into

%% Vaccuum solution

[du_zero, dv_zero, bxz_zero, byz_zero, bxy_zero, z, A_zero, d_n_zero, ...
    dummy_eig_zero,  dummy_dia_zero ,condition_A_zero] = ...
    Vacumm_plate_ij_zero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega, number_i,number_j,...
    number_m, number_n, n_thick, n_z);

%% Figures d_n vs omega

figure()
semilogy(omega/(2*pi)/omega_factor, condition_A_zero,'LineWidth',2)
xlabel('$f$ $[Hz]$','interpreter','latex','FontSize',14,'LineWidth',2)
ylabel('condition of A','interpreter','latex','FontSize',14)
hold on
grid on

%% Look for local minima and input to Secant Method and Local Minima from d_n values

dn_abs_zero = abs(d_n_zero);

TF_n_zero = islocalmin(dn_abs_zero);
TF_conzero = islocalmax(condition_A_zero);


% Obtain all local minima for possible eigenvalues

for i=1:length(omega)
    
    if (TF_n_zero(i) ==1)
        
        omega_min(i) = omega(i);
        
    end
    
end

omegas = nonzeros(omega_min);
omegas_Hz = nonzeros(omega_min)/(2*pi)/omega_factor;

%% Refine approximation of eigenfrequencies

for iter_omegas=1:length(omegas)
    
    omega0 = omegas(iter_omegas)-0.002; % 0.0001 for thin plates
    omega1 = omegas(iter_omegas)+0.002; % 0.0001 for thin plates
    %
    %     number_i= 12; % number i,j
    %     number_j= 12; % number i,j
    %
    %     number_m= 12; % number of m,n
    %     number_n= 12; % number of m,n
    
    omegas_refined = linspace(omega0,omega1, 100);
    
    [du_zero, dv_zero, bxz_zero, byz_zero, bxy_zero, z, A_zero,...
        d_n1_zero, dummy_eig1_zero,  dummy_dia1_zero,condition_A1_zero] = ...
        Vacumm_plate_ij_zero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omegas_refined, ...
        number_i,number_j, number_m, number_n, n_thick, n_z);
    
    
    TF_1_zero = islocalmin(abs(d_n1_zero));
    
    %% Obtain all local minima for possible eigenvalues
    
    for iter=1:length(omegas_refined)
        
        if (TF_1_zero(iter) ==1)
            
            omega_final(iter) = omegas_refined(iter);
            
        else
            
            omega_final(iter) =0;
            
        end
        
    end
    
    omega_final = nonzeros(omega_final)
    
    clear du_zero dv_zero dw_zero bxz_zero byz_zero bzz_zero du dv dw bzz bxz byz
    %% Mode shape reconstruction
    
    [du_zero, dv_zero, bxz_zero, byz_zero, bxy_zero, z, A_zero, d_nfinal__zero, ...
        dummy_eigfinal_zero(:,:,iter_omegas),  dummy_diafinal_zero(:,iter_omegas),condition_Afinal_zero]...
        = Vacumm_plate_ij_zero(Ef,nuf,-hf/2,hf/2,a,b,rhof, omega_final(end), ...
        number_i,number_j, number_m, number_n, n_thick, n_z);
    % NOTE: dummy_eigfinal_zero IS NOTE SORTED BASED ON THE SORTING OF EIGENVALUES
    %% analyze if we have a multiple root
    
    x= linspace(0,a,101)/hf;
    y = linspace(0,b,101)/hf;
    
    z_new = z(1,:);
    
    a_bar = a/hf;
    b_bar = b/hf;
    
    u_int_zero = nan(number_i,length(z_new));
    sigmaxz_int_zero = nan(number_i,length(z_new));
    v_int_zero = nan(number_j,length(z_new));
    sigmayz_int_zero = nan(number_j,length(z_new));
    sigmaxy_int_zero = nan(number_j,length(z_new));
    
    
    U_total_zero = nan(length(x),length(y),length(z_new));
    V_total_zero = nan(length(x),length(y),length(z_new));
    sigmaxz_total_zero = nan(length(x),length(y),length(z_new));
    sigmaxy_total_zero = nan(length(x),length(y),length(z_new));
    
    u_int2_zero = nan(number_i,length(z_new));
    sigmaxz_int2_zero = nan(number_i,length(z_new));
    v_int2_zero = nan(number_j,length(z_new));
    sigmayz_int2_zero = nan(number_j,length(z_new));
    sigmaxy_int2_zero = nan(number_j,length(z_new));
    
    
    U_total2_zero = nan(length(x),length(y),length(z_new));
    V_total2_zero = nan(length(x),length(y),length(z_new));
    sigmaxz_total2_zero = nan(length(x),length(y),length(z_new));
    sigmaxy_total2_zero = nan(length(x),length(y),length(z_new));
    
    
    eps = 10^-9;
    
    %%% i, j terms equal to zero
    
    if abs(dummy_diafinal_zero(end,iter_omegas) - dummy_diafinal_zero(end-1,iter_omegas))< eps
        
        % % NOTE: dummy_eigfinal_zero IS NOTE SORTED BASED ON THE SORTING OF EIGENVALUES

        U_root_zero = dummy_eigfinal_zero(:,end,iter_omegas); %% I changed it, put it back when done
        U_root_zero = U_root_zero/max(abs(U_root_zero));
        
        U_root2_zero = dummy_eigfinal_zero(:,end-1, iter_omegas);
        U_root2_zero = U_root2_zero/max(abs(U_root2_zero));
        
        display('multiple root')
        
        %% computation for zero terms  (i=0 or j=0); computation of displacement and stresses in the z-direction
        
        for l= 1:length(z_new)
            
            for r= 1: 2*(number_i)
                
                u_int_zero(r,l) = du_zero(r,1:2,l)*U_root_zero((r-1)*2+1:2*r);
                sigmaxz_int_zero(r,l) = bxz_zero(r,1:2,l)*U_root_zero((r-1)*2+1:2*r);
                
            end
            
        end
        
        for l= 1:length(z_new)
            
            for r= 1: 2*(number_i)
                
                v_int_zero(r,l) = dv_zero(r,1:2,l)*U_root_zero((r-1)*2+1:2*r);
                sigmayz_int_zero(r,l) = byz_zero(r,1:2,l)*U_root_zero((r-1)*2+1:2*r);
                sigmaxy_int_zero(r,l) = bxy_zero(r,1:2,l)*U_root_zero((r-1)*2+1:2*r);
                
            end
            
        end
        
        %% computation of global displacement and stresses from the contribution of all i,j
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    U_total_zero(t,o,l) = 0;
                    sigmaxz_total_zero(t,o,l) =0;
                    sigmaxy_total_zero(t,o,l) =0;
                    
                    for j=1:number_j
                        
                        r = j;
                        
                        sigmaxz_total_zero(t,o,l) = sigmaxz_total_zero(t,o,l) + sigmaxz_int_zero(r,l).*sin((j)*pi*y(o)/b_bar);
                        U_total_zero(t,o,l) = U_total_zero(t,o,l) + u_int_zero(r,l).*sin((j)*pi*y(o)/b_bar);
                        sigmaxy_total_zero(t,o,l) = sigmaxy_total_zero(t,o,l) + sigmaxy_int_zero(r,l).*cos((j)*pi*y(o)/b_bar);
  
                    end
                    
                end
                
            end
            
        end
        
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    V_total_zero(t,o,l) = 0;
                    
                    for i=1:number_i
                        
                        r = i+10;
                        
                        V_total_zero(t,o,l) = V_total_zero(t,o,l) + v_int_zero(r,l).*sin((i)*pi*x(t)/a_bar);
                        sigmaxy_total_zero(t,o,l) = sigmaxy_total_zero(t,o,l) + sigmaxy_int_zero(r,l).*cos((i)*pi*y(o)/b_bar);
                        
                    end
                    
                end
                
            end
            
        end
        
        %% the multiple case
        
        for l= 1:length(z_new)
            
            for r= 1: 2*(number_i)
                
                u_int2_zero(r,l) = du_zero(r,1:2,l)*U_root2_zero((r-1)*2+1:2*r);
                sigmaxz_int2_zero(r,l) = bxz_zero(r,1:2,l)*U_root2_zero((r-1)*2+1:2*r);
                sigmaxy_int2_zero(r,l) = bxy_zero(r,1:2,l)*U_root2_zero((r-1)*2+1:2*r);
                
            end
            
        end
        
        for l= 1:length(z_new)
            
            for r= 1: 2*(number_i)
                
                v_int2_zero(r,l) = dv_zero(r,1:2,l)*U_root2_zero((r-1)*2+1:2*r);
                sigmayz_int2_zero(r,l) = byz_zero(r,1:2,l)*U_root2_zero((r-1)*2+1:2*r);
                
            end
            
        end
        
        %% computation of global displacement and stresses from the contribution of all i,j
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    U_total2_zero(t,o,l) = 0;
                    sigmaxz_total2_zero(t,o,l) =0;
                    sigmaxy_total2_zero(t,o,l) =0;
                                            
                        for j=1:number_j
                            
                            r = j;
                            
                            sigmaxz_total2_zero(t,o,l) = sigmaxz_total2_zero(t,o,l) + sigmaxz_int2_zero(r,l).*sin((j)*pi*y(o)/b_bar);
                            U_total2_zero(t,o,l) = U_total2_zero(t,o,l) + u_int2_zero(r,l).*sin((j)*pi*y(o)/b_bar);
                            sigmaxy_total2_zero(t,o,l) = sigmaxy_total2_zero(t,o,l) + sigmaxy_int2_zero(r,l).*cos((j)*pi*y(o)/b_bar);
                            
                    end
                    
                end
                
            end
            
        end
        
        
        for l= 1:length(z_new)
            
            for t = 1:length(x)
                
                for o = 1:length(y)
                    
                    V_total2_zero(t,o,l) = 0;

                        for i=1:number_j
                            
                            r = i+10;
                            
                            V_total2_zero(t,o,l) = V_total2_zero(t,o,l) + v_int2_zero(r,l).*sin((i)*pi*x(t)/a_bar);
                            sigmaxy_total2_zero(t,o,l) = sigmaxy_total2_zero(t,o,l) + sigmaxy_int2_zero(r,l).*cos((i)*pi*y(o)/a_bar);

                    end
                    
                end
                
            end
            
        end
        
        sigmaxz_total =sigmaxz_total_zero;
        U_total = U_total_zero;
        V_total = V_total_zero;
        W_total =zeros(length(x),length(y),length(z_new));
        sigmazz_total = zeros(length(x),length(y),length(z_new));
        sigmaxx_total = zeros(length(x),length(y),length(z_new));
        W_total(1,50,1) =eps;
        sigmazz_total(1,50,1) =eps;
        sigmaxx_total(1,50,1) =eps;
        sigmaxy_total = sigmaxy_total_zero;
        
        sigmaxz_total2 = sigmaxz_total2_zero;
        U_total2 = U_total2_zero;
        V_total2 = V_total2_zero;
        W_total2 = zeros(length(x),length(y),length(z_new));
        sigmazz_total2 = zeros(length(x),length(y),length(z_new));
        sigmaxx_total2 = zeros(length(x),length(y),length(z_new));
        W_total2(1,50,1) =eps;
        sigmazz_total2(1,50,1) =eps;
        sigmaxx_total2(1,50,1) =eps;
        sigmaxy_total2 = sigmaxy_total2_zero;
        
        %% Save files
        if flag_save_dry > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            name2 = "D_intest";
            test = strcat(name2,ch2,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'W_total', 'U_total','V_total', 'sigmazz_total','sigmaxz_total', 'sigmaxx_total','sigmaxy_total', 'W_total2', 'U_total2', 'V_total2','sigmazz_total2','sigmaxz_total2', 'sigmaxx_total2', 'sigmaxy_total2','du_zero', 'dv_zero', 'bxz_zero', 'byz_zero','bxy_zero', 'U_root_zero','U_root2_zero');
            
        end
        
        if flag_save_wet > 0
            
            ch = int2str(iter_omegas);
            ch2 = int2str(a_test);
            ch3 = int2str(H_test);
            %name2 = "wetmodes";
            name2 = "W_in";
            test = strcat(name2,ch2,ch3,ch,".mat");
            
            save(test,'x','y', 'z_new','omega_final', 'W_total', 'U_total', 'V_total','sigmazz_total','sigmaxz_total', 'W_total2', 'U_total2','V_total2', 'sigmazz_total2','sigmaxz_total2', 'du_zero', 'dv_zero', 'bxz_zero', 'byz_zero', 'U_root_zero','U_root2_zero');
            
        end
        
    end
    
    clear omega_final
end