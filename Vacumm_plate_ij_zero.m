%% Function to compute material constants, A matrix, approximate eigenfrequencies

function [du, dv, bxz, byz, bxy, z, Final_A, d_n, dummy_eig, dummy_dia, condition_A] = ...
    Vacumm_plate_ij_zero_v2(E,nu,hinit,hfinal,a,b,rho, omega, number_i, number_j, number_m, number_n, n_thick, n_z)

h = hfinal - hinit;
a_bar = a/h;
b_bar = b/h;

%% Computation of trigonometric integrals

for Im = 1:number_m
    
    for Ii = 1:number_i
        
        %%%sin(x)*sin(x)
        
        if (Ii == 0) || (Im == 0) % THIS CONDITION IS NEVER TRUE
            
            Ix(Ii,Im) = 0;
            
        else
            
            if   (Ii == Im)
                
                Ix(Ii,Im) = a_bar/2;
                
            else
                
                
                Ix(Ii,Im) = 0;
                
            end
        end
        
        
    end
end

for In = 1:number_n
    
    for Ij = 1:number_j
        
        
        %%%sin(y)*sin(y)
        
        if (Ij == 0) || (In == 0)
            
            Iy(Ij,In) = 0;
            
        else
            
            if   (Ij == In)
                
                Iy(Ij,In) = b_bar/2;
                
            else
                
                
                Iy(Ij,In) = 0;
                
            end
        end
        
    end
    
end

%%%1*cos(x)
for Im = 1:1 % WHY ADDING A FOR LOOP?
    
    for Ii=1:1
        
        if (Ii-1 == 0) && (Im-1 == 0)
            
            Icosx(Ii,Im) = a_bar;
            
        else
            
            if   (Ii-1 == Im-1)
                
                Icosx(Ii,Im) = a_bar/2;
                
            else
                
                
                Icosx(Ii,Im) = 0;
                
            end
        end
        
        
    end
end

for In = 1:1
    
    for Ij=1:1
        
        %%%cos(y)*cos(y)
        
        if (Ij-1 == 0) &&  (In-1 ==0)
            
            Icosy(Ij,In) = b_bar;
            
        else
            
            if   (Ij-1 == In-1)
                
                Icosy(Ij,In) = b_bar/2;
                
            else
                
                Icosy(Ij,In) = 0;
                
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

for Iomega=1:length(omega)
    
%     it=Iomega
    
    %% case of i =0
    for Ii = 1:1
        
        for Ij= 1:number_j
            
            r = (Ii-1)*number_i+Ij;
            
            [du(r,:,:), ~, ~, ~, bxz(r,:,:), ~,~,~,bxy1(r,:,:), z(r,:,:)] = new_pagano_non(Ef,nuf,Ii-1, Ij, hinit, hfinal,a,b,rhof, omega(Iomega));
            
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
    
    for Ii = 1:number_i
        
        for Ij=  1:1
            
            r = (Ij-1)*number_i+Ii;
            
            [~, dv(r,:,:), ~, ~, ~, byz(r,:,:),~,~,bxy2(r,:,:), z(r,:,:)] = new_pagano_non(Ef,nuf,Ii, Ij-1, hinit, hfinal,a,b,rhof, omega(Iomega));
            
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
    
    Final_A1 = [eq2i'
        eq5i'];
    
    Final_A2 = [eq3j'
        eq6j'];

    Final_A = [Final_A1 zeros(2*number_i,2*number_j)
        zeros(2*number_i,2*number_j) Final_A2];
    % assembly of boundary conditions
    
    %% Compute eigenvalues of the A matrix

    [~,dummy,dummy_eig] = svd(Final_A);
    
    dummy_dia = diag(dummy);

    [dummy_dia, ind] = sort(abs(dummy_dia),'descend');
    
    %% get eigenvectors

     dummy_eig = dummy_eig(:,ind); %sort the eigenvectors
    
    condition_A(Iomega) = cond(Final_A); % condition number of the whole matrix with u and v components
    d_n(Iomega) = dummy_dia(end); % smallest eigenvalue of A(omega)
    
    
    %% through-the-thickness basis functions
    
    du = [du 
        zeros(number_i,6,n_thick)];
    
    dv = [zeros(number_i,6,n_thick) 
        dv];
    
    bxz = [bxz 
        zeros(number_i,6,n_thick)];
    
    byz = [zeros(number_i,6,n_thick) 
        byz];
    
    bxy = [bxy1 
        bxy2];
    
    
end

