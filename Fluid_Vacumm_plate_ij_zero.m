%% Function to compute material constants, A matrix, approximate eigenfrequencies

function[du, dv, bxz, byz, bxy, z, Final_A, d_n, dummy_eig, dummy_dia, condition_A] = Fluid_Vacumm_plate_ij_zero(E,nu,hinit,hfinal,a,b,rho, omega, number_i, number_j, number_m, number_n, n_thick, n_z)

h = hfinal - hinit;
a_bar = a/h;
b_bar = b/h;

%% Computation of trigonometric integrals

for m = 1:number_m
    
    for i = 1:number_i
        
        %%%sin(x)*sin(x)
        
        if (i == 0) || (m == 0)
            
            Ix(i,m) = 0;
            
        else
            
            if   (i == m)
                
                Ix(i,m) = a_bar/2;
                
            else
                
                
                Ix(i,m) = 0;
                
            end
        end
        
        
    end
end
        
for n = 1:number_n
    
    for j = 1:number_j
        
        
        %%%sin(y)*sin(y)
        
        if (j == 0) || (n == 0)
            
            Iy(j,n) = 0;
            
        else
            
            if   (j == n)
                
                Iy(j,n) = b_bar/2;
                
            else
                
                
                Iy(j,n) = 0;
                
            end
        end
        
    end
    
end

%%%1*cos(x)
for m = 1:1
    
    for i=1:1

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

for n = 1:1
    
    for j=1:1
        
        %%%cos(y)*cos(y)
        
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

%% Computation of Krs and Mrs matrices

                K_xz= Icosx*Iy;
                K_yz = Ix*Icosy;

%% computation of W-ij terms and sigma-ij terms for boundary conditions
Ef= E;
nuf = nu;
rhof = rho;

bxz = zeros(number_i,6,n_thick);
byz = zeros(number_i,6,n_thick);
bxy = zeros(number_i,6,n_thick);
z = zeros(number_i,1,n_thick);

sigmaxz_ij_top = zeros(number_i,6);
sigmayz_ij_top = zeros(number_i,6);
sigmaxz_ij_bottom = zeros(number_i,6);
sigmayz_ij_bottom = zeros(number_i,6);

d_n = nan(length(omega),1);
condition_A = nan(length(omega),1);

for t=1:length(omega)
    
    it=t
    
    %% case of i =0
    for i = 1:1
        
        for j= 1:number_j
            
            r = (i-1)*number_i+j;
            
            [du(r,:,:), ~, ~, ~, bxz(r,:,:), ~, ~,~,bxy(r,:,:), z(r,:,:)] = new_pagano_non(Ef,nuf,i-1, j, hinit, hfinal,a,b,rhof, omega(t));

            sigmaxz_ij_top(r,:) = bxz(r,:,n_thick);
            sigmaxz_ij_bottom(r,:) = bxz(r,:,1);

        end
        
    end
    
    %% Computation of \bar{M} and \bar{K}
    
    for s = 1:number_m
        
        for r=1:number_i
            
            for l = 1:2 %n_z %number of unknowns through thickness dimension W(z), sigma(z)
                
                %%% for eq2
                
                K_final2((r-1)*2+l,s)= K_xz(r,s)*sigmaxz_ij_top(r,l);  %sigma_xz = 0

                %%% for eq5
                
                K_final5((r-1)*2+l,s) = K_xz(r,s)*sigmaxz_ij_bottom(r,l);  %sigma_xz = 0
                
            end
        end
        
    end
    
    
    eq2i = K_final2;
    eq5i= K_final5;
    
    %% case of j=0
    
    for i = 1:number_i
        
        for j=  1:1
                
                r = (j-1)*number_i+i;
                
                [~, dv(r,:,:), ~, ~, ~, byz(r,:,:),~,~,bxy(r,:,:), z(r,:,:)] = new_pagano_non(Ef,nuf,i, j-1, hinit, hfinal,a,b,rhof, omega(t));
             
                sigmayz_ij_top(r,:) = byz(r,:,n_thick);
                sigmayz_ij_bottom(r,:) = byz(r,:,1);

        end
        
    end
    
    %% Computation of \bar{M} and \bar{K}
    
    for s = 1:number_m
        
        for r=1:number_i
            
            for l = 1:2 %n_z %number of unknowns through thickness dimension W(z), sigma(z)
                
                %%% for eq3
                
                K_final3((r-1)*2+l,s) = K_yz(r,s)*sigmayz_ij_top(r,l);  %sigma_yz = 0
                
                %%% for eq6
                
                K_final6((r-1)*2+l,s)= K_yz(r,s)*sigmayz_ij_bottom(r,l);  %sigma_yz = 0
                
            end
        end
        
    end
    
    
    
    eq3j = K_final3;
    eq6j= K_final6;
    
    %% Assembly of i=0 and j=0 cases into the stiffness matrix
    
    Final_A = [eq2i'
        eq5i'
        eq3j'
        eq6j']; % assembly of boundary conditions
    
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