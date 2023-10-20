%% Function to compute material constants, A matrix, approximate eigenfrequencies

function[du, dv, dw, bzz, bxz, byz, bxx, byy, bxy, z, Final_A, d_n, dummy_eig, dummy_dia, condition_A] = Fluid_Vacumm_plate_nonzero(E,nu,hinit,hfinal,a,b,rho, omega, rho_fluid, H, number_i, number_j, number_m, number_n, n_thick, n_z)

I1 = nan(number_i,number_m);
I2 = nan(number_i,number_m);
I1_2 = nan(number_j,number_n);
I3 = nan(number_j,number_n);

h = hfinal - hinit;
a_bar = a/h;
b_bar = b/h;

H_bar = H/h;

% % % % % 
% % % % % %% Computation of trigonometric integrals
% % % % % 
% % % % % for m = 1:number_m
% % % % %     
% % % % %     for n=1:number_n
% % % % %         
% % % % %         for i = 1:number_i
% % % % %             
% % % % %             for j=1:number_j
% % % % %                 
% % % % %                 %%%sin(x)*cos(x)
% % % % %                 
% % % % %                 if (mod(i,2) + mod((m-1),2)) ==1
% % % % %                     
% % % % %                     I1((i-1)*number_i+j,(m-1)*number_m+n) = a_bar*(1/(i+(m-1))+1/(i-(m-1)))/(pi);
% % % % %                     
% % % % %                 else
% % % % %                     
% % % % %                     I1((i-1)*number_i+j,(m-1)*number_m+n) = 0;
% % % % %                     
% % % % %                 end
% % % % %                 
% % % % %                 %%%sin(y)*cos(y)
% % % % %                 
% % % % %                 if  (mod(j,2) + mod((n-1),2)) == 1
% % % % %                     
% % % % %                     I1_2((i-1)*number_i+j,(m-1)*number_m+n) = b_bar*(1/(j+(n-1))+1/(j-(n-1)))/(pi);
% % % % %                     
% % % % %                 else
% % % % %                     
% % % % %                     I1_2((i-1)*number_i+j,(m-1)*number_m+n) = 0;
% % % % %                     
% % % % %                 end
% % % % %                 
% % % % %                 %%%cos(x)*cos(x)
% % % % %                 
% % % % %                 if (i == (m))
% % % % %                     
% % % % %                     I2((i-1)*number_i+j,(m-1)*number_m+n) = a_bar/2;
% % % % %                 else
% % % % %                     
% % % % %                     I2((i-1)*number_i+j,(m-1)*number_m+n) = 0;
% % % % %                     
% % % % %                 end
% % % % %                 
% % % % %                 %%%cos(y)*cos(y)
% % % % %                 
% % % % %                 if (j == (n))
% % % % %                     
% % % % %                     I3((i-1)*number_i+j,(m-1)*number_m+n) = b_bar/2;
% % % % %                     
% % % % %                 else
% % % % %                     
% % % % %                     I3((i-1)*number_i+j,(m-1)*number_m+n) = 0;
% % % % %                     
% % % % %                 end
% % % % %                 
% % % % %             end
% % % % %             
% % % % %         end
% % % % %         
% % % % %         %% Computation of epsilon_mn
% % % % %         
% % % % %         if (m-1) == 0 && (n-1) ==0
% % % % %             
% % % % %             epsilon((m-1)*number_m+n)=1;
% % % % %             
% % % % %             
% % % % %         elseif ((m-1)~=0 && (n-1)==0)
% % % % %             
% % % % %             epsilon((m-1)*number_m+n) =2;
% % % % %             
% % % % %         elseif ((n-1)~=0 && (m-1)==0)
% % % % %             
% % % % %             epsilon((m-1)*number_m+n) =2;
% % % % %             
% % % % %         else
% % % % %             
% % % % %             epsilon((m-1)*number_m+n) =4;
% % % % %             
% % % % %         end
% % % % %         
% % % % %         %% Computation of zeta_mn
% % % % %         
% % % % %         if (m-1) ==0 && (n-1) ==0
% % % % %             
% % % % %             zeta((m-1)*number_m+n) = 1;
% % % % %             
% % % % %         elseif ((m-1) ~=0 && (n-1) ==0)
% % % % %             
% % % % %             zeta((m-1)*number_m+n) =  1/2;
% % % % %             
% % % % %         elseif  ((n-1) ~=0 && (m-1) ==0)
% % % % %             
% % % % %             zeta((m-1)*number_m+n) =  1/2;
% % % % %             
% % % % %         else
% % % % %             
% % % % %             zeta((m-1)*number_m+n) = 1/4;
% % % % %             
% % % % %             
% % % % %         end
% % % % %         
% % % % %         %% Computation of gamma_mn, A_mn and C_mn
% % % % %         
% % % % %         gamma((m-1)*number_m+n) = sqrt(((m-1)*pi/a_bar)^2+((n-1)*pi/b_bar)^2);
% % % % %         
% % % % %         if (m-1) == 0 && (n-1) == 0
% % % % %             
% % % % %             Amn((m-1)*number_m+n) = epsilon((m-1)*number_m+n)/(a_bar*b_bar);
% % % % %             
% % % % %             Cmn((m-1)*number_m+n) = - rho_fluid/rho*H*Amn((m-1)*number_m+n);
% % % % %             
% % % % %         else
% % % % %             
% % % % %             Amn((m-1)*number_m+n) = epsilon((m-1)*number_m+n)/(a_bar*b_bar*gamma((m-1)*number_m+n));
% % % % %             
% % % % %             Cmn((m-1)*number_m+n) = - rho_fluid/rho*tanh(gamma((m-1)*number_m+n)*H)*Amn((m-1)*number_m+n);
% % % % %             
% % % % %         end
% % % % %         
% % % % %         
% % % % %     end
% % % % %     
% % % % % end
% % % % % 
% % % % % %% Computation of Krs and Mrs matrices
% % % % % 
% % % % % K_zz = I1.*I1_2;
% % % % % K_xz = I2.*I1_2;
% % % % % K_yz = I1.*I3;
% % % % % 
% % % % % for m = 1:number_m
% % % % %     
% % % % %     for n=1:number_n
% % % % %         
% % % % %         for i = 1:number_i
% % % % %             
% % % % %             for j=1:number_j
% % % % %                 
% % % % %                 M((i-1)*number_i+j,(m-1)*number_m+n) =  Cmn((m-1)*number_m+n)*zeta((m-1)*number_m+n)*a_bar*b_bar*K_zz((i-1)*number_i+j,(m-1)*number_m+n);
% % % % %                 
% % % % %             end
% % % % %             
% % % % %         end
% % % % %         
% % % % %     end
% % % % %     
% % % % % end
% % % % % 






%% Computation of trigonometric integrals

for m = 1:number_m
    
    for i = 1:number_i
        
        %%%sin(x)*sin(x)
        
        if (i == 0) || (m == 0)
            
            I1(i,m) = 0;
            
        else
            
            if   (i == m)
                
                I1(i,m) = a_bar/2;
                
            else
                
                
                I1(i,m) = 0;
                
            end
        end
        
        
        %%%cos(x)*cos(x)
        
        if (i == 0) && (m == 0)
            
            I2(i,m) = a_bar;
            
        else
            
            if   (i == m)
                
                I2(i,m) = a_bar/2;
                
            else
                
                
                I2(i,m) = 0;
                
            end
        end
        
        
    end
end


for n = 1:number_n
    
    for j = 1:number_j
        
        
        %%%sin(y)*sin(y)
       
        if (j == 0) || (n == 0)
            
            I1_2(j,n) = 0;
            
        else
            
            if   (j == n)
                
                I1_2(j,n) = b_bar/2;
                
            else
                
                
                I1_2(j,n) = 0;
                
            end
        end
        
        
        
        %%%cos(y)*cos(y)
        
        if (j == 0) &&  (n ==0)
            
            I3(j,n) = b_bar;
            
        else
            
            if   (j == n)
                
                I3(j,n) = b_bar/2;
                
            else
                
                I3(j,n) = 0;
                
            end
            
        end
        
    end
    
end

%% integrals for Imn term (sin*cos)

Ix = nan(number_i,number_m);
Iy = nan(number_j,number_n);

for m = 1:number_m
    
    for i = 1:number_i
        
        %%%sin(x)*cos(x)
        
        if  (i) == 0
            
            Ix(i,m) = 0;
            
        else
            
            if (mod((i),2) + mod((m-1),2)) == 1
                
                Ix(i,m) = a_bar*(1/((i)+(m-1))+1/((i)-(m-1)))/(pi);
                
            else
                
                Ix(i,m) = 0;
                
            end
            
        end
        
        
    end
    
end



for n = 1:number_n
    
    for j = 1:number_j
        
        %%%sin(y)*cos(y)
        
        if  (j) == 0
            
            Iy(j,n) = 0;
            
        else
            
            
            if  (mod((j),2) + mod((n-1),2)) == 1
                
                Iy(j,n) = b_bar*(1/((j)+(n-1))+1/((j)-(n-1)))/(pi);
                
            else
                
                Iy(j,n)= 0;
                
            end
            
        end
        
    end
    
end


%% %% integrals for Imn term (sin*cos) version 2

Ix2 = nan(number_i,number_m);
Iy2 = nan(number_j,number_n);

for m = 1:number_m
    
    for i = 1:number_i
        
        %%%sin(x)*cos(x)
        
        if  (i) == 0
            
            Ix2(i,m) = 0;
            
        else
            
            if (mod((i),2) + mod((m),2)) == 1
                
                Ix2(i,m) = a_bar*(1/((i)+(m))+1/((i)-(m)))/(pi);
                
            else
                
                Ix2(i,m) = 0;
                
            end
            
        end
        
        
    end
    
end



for n = 1:number_n
    
    for j = 1:number_j
        
        %%%sin(y)*cos(y)
        
        if  (j) == 0
            
            Iy2(j,n) = 0;
            
        else
            
            
            if  (mod((j),2) + mod((n),2)) == 1
                
                Iy2(j,n) = b_bar*(1/((j)+(n))+1/((j)-(n)))/(pi);
                
            else
                
                Iy2(j,n)= 0;
                
            end
            
        end
        
    end
    
end

%%%%%



%% integrals for Icoscos version 2

for m = 1:number_m
    
    for i = 1:number_i
        
        %%%cos(x)*cos(x)
        
        if (i-1 == 0) && (m-1 == 0)
            
            Icosx(i,m) = a_bar;
            
        else
            
            if   (i-1 == m-1)
                
                Icosx(i,m) = a_bar/2;
                
            else
                
                
                Icosx(i,m) = 0;
                
            end
        end
        
    end
    
end

        %%%cos(y)*cos(y)
  
        
for n = 1:number_n
    
    for j = 1:number_j
        
        if (j-1 == 0) &&  (n-1 ==0)
            
            Icosy(j,n) = b_bar;
            
        else
            
            if   (j-1 == n-1)
                
                Icosy(j,n) = b_bar/2;
                
            else
                
                Icosy(j,n) = 0;
                
            end
            
        end
        
    end
    
end
        
        %%%%%%

gamma = nan(number_m*number_n,1);
Cmn = nan(number_m*number_n,1);

for m = 1:number_m
    
    for n = 1:number_n   
        
        %% Computation of gamma_mn, A_mn and C_mn
        
        s = (m-1)*number_m+n;
        
        gamma(s) = sqrt(((m-1)*pi/a_bar)^2+((n-1)*pi/b_bar)^2);
        
        if (m-1) == 0 && (n-1) == 0
            
            Cmn(s) =  -rho_fluid/rho*(-H_bar);
            
        else

            Cmn(s) =  -rho_fluid/rho*(-tanh(gamma(s)*H_bar))/gamma(s);
            
        end
        
        
    end
    
end

%% Computation of Krs and Mrs matrices

K_zz = nan(number_i*number_j,number_m*number_n);
K_xz = nan(number_i*number_j,number_m*number_n);
K_yz = nan(number_i*number_j,number_m*number_n);
Imn = nan(number_i*number_j,number_m*number_n);
Imn2 = nan(number_i*number_j,number_m*number_n);
Icoscos = nan(number_i*number_j,number_m*number_n);

for m = 1:number_m
    
    for n = 1:number_n
        
        for i = 1:number_i
            
            for j = 1:number_j
                
                r = (i-1)*number_i+j;
                s = (m-1)*number_m+n;

                K_xz(r,s)= I2(i,m)*I1_2(j,n); % both start from 1
                K_yz(r,s) = I1(i,m)*I3(j,n); % both start from 1
                K_zz(r,s) = I1(i,m)*I1_2(j,n); % both start from 1

                
                Imn(r,s) = Ix(i,m)*Iy(j,n); % for sigmazz_top i,j start in 1 but m,n in 0
                Imn2(r,s) = Ix2(i,m)*Iy2(j,n); %i,j and m,n starting both from 1 (Amn term)
                Icoscos(r,s) = Icosx(i,m)*Icosy(j,n); % here all index should start from zero
                
            end
            
        end
    end
end

%% computation of W-ij terms and sigma-ij terms for boundary conditions
Ef= E;
nuf = nu;
rhof = rho;

du = nan(number_i*number_j,6,n_thick);
dv = nan(number_i*number_j,6,n_thick);
dw = nan(number_i*number_j,6,n_thick);
bzz = nan(number_i*number_j,6,n_thick);
bxz = nan(number_i*number_j,6,n_thick);
byz = nan(number_i*number_j,6,n_thick);
bxx = nan(number_i*number_j,6,n_thick);
byy = nan(number_i*number_j,6,n_thick);
bxy = nan(number_i*number_j,6,n_thick);
z = nan(number_i*number_j,1,n_thick);

dw_ij_top = nan(number_i*number_j,6);
sigmazz_ij_top = nan(number_i*number_j,6);
sigmaxz_ij_top = nan(number_i*number_j,6);
sigmayz_ij_top = nan(number_i*number_j,6);
sigmazz_ij_bottom = nan(number_i*number_j,6);
sigmaxz_ij_bottom = nan(number_i*number_j,6);
sigmayz_ij_bottom = nan(number_i*number_j,6);

d_n = nan(length(omega),1);
condition_A = nan(length(omega),1);

for t=1:length(omega)
    
    it=t
    for i = 1:number_i
        
        for j= 1:number_j
            
            r = (i-1)*number_i+j;
            
            [du(r,:,:), dv(r,:,:), dw(r,:,:), bzz(r,:,:), bxz(r,:,:), byz(r,:,:), bxx(r,:,:),byy(r,:,:),bxy(r,:,:), z(r,:,:)] = new_pagano_non(Ef,nuf,i, j, hinit, hfinal,a,b,rhof, omega(t));
            
            dw_ij_top(r,:) = dw(r,:,n_thick);
            sigmazz_ij_top(r,:) = bzz(r,:,n_thick);
            sigmaxz_ij_top(r,:) = bxz(r,:,n_thick);
            sigmayz_ij_top(r,:) = byz(r,:,n_thick);
            
            sigmazz_ij_bottom(r,:) = bzz(r,:,1);
            sigmaxz_ij_bottom(r,:) = bxz(r,:,1);
            sigmayz_ij_bottom(r,:) = byz(r,:,1);
            
        end
        
    end
    
    %% Computation of \bar{M} and \bar{K}
    
    M_final = nan((number_i)*number_j*6,number_m*number_n);
    K_final = nan((number_i)*(number_j)*6,number_m*number_n);
    K_final2 = nan((number_i)*(number_j)*6,number_m*number_n);
    K_final3 = nan((number_i)*(number_j)*6,number_m*number_n);
    K_final4 = nan((number_i)*(number_j)*6,number_m*number_n);
    K_final5 = nan((number_i)*(number_j)*6,number_m*number_n);
    K_final6 = nan((number_i)*(number_j)*6,number_m*number_n);
    
    M = nan(number_i*number_j,number_m*number_n);

for m = 1:number_m
    
    for n=1:number_n
        
        for i = 1:number_i
            
            for j=1:number_j
                
                r = (i-1)*number_i+j;
                s = (m-1)*number_m+n;
                
                
                    M(r,s) = Cmn(s)*Imn(r,s);

                
            end
            
        end
        
    end
    
end
    
    
    for s = 1:number_m*number_n
        
        for r=1:number_i*number_j
            
            for l = 1:n_z %number of unknowns through thickness dimension W(z), sigma(z)
                
                M_final((r-1)*n_z+l,s) = M(r,s)*dw_ij_top(r,l);  %sigma_zz = p
                
                K_final((r-1)*n_z+l,s) = Imn(r,s)*sigmazz_ij_top(r,l);  %sigma_zz = p
                
                %%% for eq2
                
                K_final2((r-1)*n_z+l,s)= K_xz(r,s)*sigmaxz_ij_top(r,l);  %sigma_xz = 0
                
                %%% for eq3
                
                K_final3((r-1)*n_z+l,s) = K_yz(r,s)*sigmayz_ij_top(r,l);  %sigma_yz = 0
                
                %%% for eq4 % K_zz(r,s)
                
                K_final4((r-1)*n_z+l,s)= K_zz(r,s)*sigmazz_ij_bottom(r,l); %sigma_zz = 0
                
                %%% for eq5
                
                K_final5((r-1)*n_z+l,s) = K_xz(r,s)*sigmaxz_ij_bottom(r,l);  %sigma_xz = 0
                
                %%% for eq6
                
                K_final6((r-1)*n_z+l,s)= K_yz(r,s)*sigmayz_ij_bottom(r,l);  %sigma_yz = 0
                
            end
        end
        
    end
    
    
    %% Eigenvalue problem, matrix assembly [from i,j=(1,1)]
    
    eq1 = K_final - omega(t)^2*M_final; % impose the 6 boundary condition for plate
    eq2 = K_final2;
    eq3 = K_final3;
    eq4= K_final4;
    eq5= K_final5;
    eq6= K_final6;
    
    Final_A = [eq1'
        eq2'
        eq3'
        eq4'
        eq5'
        eq6']; % assembly of boundary conditions

    %% Compute eigenvalues of the A matrix
    
    [~,dummy,dummy_eig] = svd(Final_A);
    
    dummy_dia = diag(dummy);
    [dummy_dia, ind] = sort(abs(dummy_dia),'descend');
    dummy_eig = dummy_eig(:,ind);
    
    condition_A(t) = cond(Final_A);
    d_n(t) = dummy_dia(end); % smallest eigenvalue of A(omega)
    %d_n(t) = eigs(Final_A,1,'smallestabs');
    
end

end