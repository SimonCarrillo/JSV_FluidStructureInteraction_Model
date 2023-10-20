clc
clear all

%% freq from COMSOl

hf = 0.02;

E= 69*10^9; % in Pa
nu= 0.3;
rho = 2.702*10^3; % in kg/m3

rhof = 1000;

tilde_rho = rhof/rho;

h = 0.02;

omega_factor = sqrt(rho/E)*h;

flag_save = 1;

%% read comsol freq

comsol_f1 = load('Param_a2_comsol.txt');
comsol_f2 = load('Param_a5_comsol.txt');
comsol_f3 = load('Param_a10_comsol.txt');
comsol_f4 = load('Param_a20_comsol.txt');


for i =1:4
    
    comsol_non1(:,i) = comsol_f1((i-1)*15+1:i*15,3);
    comsol_non2(:,i) = comsol_f2((i-1)*15+1:i*15,3);
    comsol_non3(:,i) = comsol_f3((i-1)*15+1:i*15,3);
    comsol_non4(:,i) = comsol_f4((i-1)*15+1:i*15,3);
    
end

eps = 0.0001;
for i= 1:4
    for t = 1:14

        if abs(comsol_non1(t,i) - comsol_non1(t+1,i)) < eps
            
            comsol_non1(t,i) = 0;
            
        end

    end
    comsol_1(:,i) = nonzeros(comsol_non1(:,i));
end


for i= 1:4
    for t = 1:14

        if abs(comsol_non2(t,i) - comsol_non2(t+1,i)) < eps
            
            comsol_non2(t,i) = 0;
            
        end

    end
    comsol_2(:,i) = nonzeros(comsol_non2(:,i));
end

for i= 1:4
    for t = 1:14

        if abs(comsol_non3(t,i) - comsol_non3(t+1,i)) < eps
            
            comsol_non3(t,i) = 0;
            
        end

    end
    comsol_3(:,i) = nonzeros(comsol_non3(:,i));
end


for i= 1:4
    for t = 1:14

        if abs(comsol_non4(t,i) - comsol_non4(t+1,i)) < eps
            
            comsol_non4(t,i) = 0;
            
        end

    end
    comsol_4{:,i} = nonzeros(comsol_non4(:,i));
end


%% read theory results FSI

for a_test = 1:4
    
    for H_test = 1:4
        
        ch2 = int2str(a_test);
        ch3 = int2str(H_test);
        name2 = "Parametric_";
        test = strcat(name2,ch2,ch3,".mat");
        
        
        theory{a_test,H_test,:}= load(test).omegas;
        
    end
end

%% error computation

for H_test = 1:4
    
    error1(:,H_test) = abs((comsol_1(1:9,H_test) - theory{1,H_test}(1:9))./comsol_1(1:9,H_test))*100
    error2(:,H_test) = abs((comsol_2(1:9,H_test) - theory{2,H_test}(1:9))./comsol_2(1:9,H_test))*100
    error3(:,H_test) = abs((comsol_3(1:9,H_test) - theory{3,H_test}(1:9))./comsol_3(1:9,H_test))*100
    error4(:,H_test) = abs((comsol_4{H_test}(1:8) - theory{4,H_test}(1:8))./comsol_4{H_test}(1:8))*100

end


%% NAVMI factors

%% read theory results vacuum

for a_test = 1:4
    
    for H_test = 1:1
        
        ch2 = int2str(a_test);
        ch3 = int2str(H_test);
        name2 = "Parametricdry_";
        test = strcat(name2,ch2,ch3,".mat");
        
        
        theory_vacuum{a_test,:}= load(test).omegas;
        
    end
end

a_parametric = [0.04 0.1 0.2 0.4];
H_parametric = [0.04 0.1 0.2 0.4];


for a_test = 1:4
    
    tilde_a = a_parametric(a_test)/h;
    
    for H_test = 1:4

        NAVMI{a_test,H_test,:} = ((theory_vacuum{a_test,:}(1:8)./theory{a_test,H_test,:}(1:8)).^2-1)./(tilde_rho*tilde_a);

    end

end



for a_test = 1:1
    
     for H_test = 1:1
         

H(H_test,:) = NAVMI{a_test,H_test,:}(1:8);
H(H_test+1,:) = NAVMI{a_test,H_test+1,:}(1:8);
H(H_test+2,:) = NAVMI{a_test,H_test+2,:}(1:8);
H(H_test+3,:) = NAVMI{a_test,H_test+3,:}(1:8);

H2(H_test,:) = NAVMI{a_test+1,H_test,:}(1:8);
H2(H_test+1,:) = NAVMI{a_test+1,H_test+1,:}(1:8);
H2(H_test+2,:) = NAVMI{a_test+1,H_test+2,:}(1:8);
H2(H_test+3,:) = NAVMI{a_test+1,H_test+3,:}(1:8);

H3(H_test,:) = NAVMI{a_test+2,H_test,:}(1:8);
H3(H_test+1,:) = NAVMI{a_test+2,H_test+1,:}(1:8);
H3(H_test+2,:) = NAVMI{a_test+2,H_test+2,:}(1:8);
H3(H_test+3,:) = NAVMI{a_test+2,H_test+3,:}(1:8);

H4(H_test,:) = NAVMI{a_test+3,H_test,:}(1:8);
H4(H_test+1,:) = NAVMI{a_test+3,H_test+1,:}(1:8);
H4(H_test+2,:) = NAVMI{a_test+3,H_test+2,:}(1:8);
H4(H_test+3,:) = NAVMI{a_test+3,H_test+3,:}(1:8);

     end
     
end


for H_test = 2:2
    
     for a_test = 1:1
         

a(a_test,:) = NAVMI{a_test,H_test,:}(1:8);
a(a_test+1,:) = NAVMI{a_test+1,H_test,:}(1:8);
a(a_test+2,:) = NAVMI{a_test+2,H_test,:}(1:8);
a(a_test+3,:) = NAVMI{a_test+3,H_test,:}(1:8);

     end
     
end


for i =1:2
    
    figure(10+i)
    leg = plot(H_parametric/h,H(:,i),'o-','LineWidth',4)
     hold on
    leg2 = plot(H_parametric/h,H2(:,i),'o-','LineWidth',4)
    leg3 = plot(H_parametric/h,H3(:,i),'o-','LineWidth',4)
    leg4 = plot(H_parametric/h,H4(:,i),'o-','LineWidth',4)
     box on
    ax = gca
    ax.LineWidth = 2
    set(gca,'fontname','times');
    set(gca,'Fontsize',30);
    axis tight
    ax.FontSize = 30;
    
    if i==1
    legend([leg leg2 leg3 leg4],'$\tilde{a}=2$','$\tilde{a}=5$','$\tilde{a}=10$','$\tilde{a}=20$','interpreter','latex', 'Location','northwest')
    end
    xlabel('$\tilde{H}$','interpreter','latex')
    ylabel('$\Gamma$','interpreter','latex')
    hold off
    
        
    if flag_save > 0
    
    fighdl = figure(10+i);
    set(fighdl,'Units','Inches');
    pos = get(fighdl,'Position');
    set(fighdl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(fighdl, ['Fig_namvi',num2str(i)],'-dpdf','-r0');
    %exportgraphics(fighdl,['Fig_2a',ch],'.pdf', );
    
    end
%     
   
end
 box on
    ax = gca
    ax.LineWidth = 2
    set(gca,'fontname','times');
    set(gca,'Fontsize',30);
    axis tight
    ax.FontSize = 30;
    legend([leg leg2 leg3 leg4],'$\tilde{a}=2$','$\tilde{a}=5$','$\tilde{a}=10$','$\tilde{a}=20$','interpreter','latex', 'Location','best')
    xlabel('$\tilde{H}$','interpreter','latex')
    ylabel('$\Gamma$','interpreter','latex')
    hold off
%     
% figure()
% for i =1:4
%     
%     leg(i) = plot(a_parametric/h,a(:,i),'o-','LineWidth',4)
%     hold on
%     
% end
%  box on
%     ax = gca
%     ax.LineWidth = 2
%     set(gca,'fontname','times');
%     set(gca,'Fontsize',30);
%     ax.FontSize = 30;
%     legend([leg(1) leg(2) leg(3) leg(4)],'$1^{st}$','$2^{nd}$','$3^{rd}$', '$4^{th}$','interpreter','latex', 'Location','northwest')
%     xlabel('$\tilde{a}$','interpreter','latex')
%     ylabel('NAVMI','interpreter','latex')
%     hold off