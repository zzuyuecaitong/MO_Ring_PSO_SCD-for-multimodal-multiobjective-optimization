load data_of_fig5 
% generate figure
figure_FontSize=15;
Mean_DP_mo_ring_pso_scd=mean(DP_mo_ring_pso_scd,3);
Mean_DP_Omni_Opt=mean(DP_Omni_Opt,3);
Mean_DP_Decision_niched_NSGAII=mean(DP_Decision_niched_NSGAII,3);
figure
plot(1:100,Mean_DP_mo_ring_pso_scd(:,1),'b-','LineWidth',2)
hold on
plot(1:100,Mean_DP_mo_ring_pso_scd(:,2),'g:','LineWidth',2)
plot(1:100,Mean_DP_mo_ring_pso_scd(:,3),'r--','LineWidth',2)
plot(1:100,Mean_DP_mo_ring_pso_scd(:,4),'k-.','LineWidth',2);
hold off
axis([0 100 0.23 0.28]);

set(gca,'Fontname','Times New Roman','Fontsize',figure_FontSize)

P1=legend ('Region 1', 'Region 2' ,'Region 3' ,'Region 4');

 L_xoff   =  0.64    ;
 L_yoff   =  0.69    ;
 L_width  =  0.12  ;
 L_height =  0.001  ;
    
        set(P1,'Position',[L_xoff,L_yoff,L_width, L_height])
        set(P1,'FontName','Times New roman','FontSize',figure_FontSize-5 )

title ('MO-Ring-PSO-SCD','FontName','Times New roman','Fontsize',figure_FontSize)
xlabel ('Generation number','FontName','Times New roman','Fontsize',figure_FontSize)
ylabel('Proportion','FontName','Times New roman','Fontsize',figure_FontSize)

set(gcf,'Units','centimeters');%把图像大小单位用厘米衡量

set(gcf,'Position',[5 5  12.9 10.5]);
set(gca,'Position',[.15 .12 .81 .79]);




figure
plot(1:100,Mean_DP_Omni_Opt(:,1),'b-','LineWidth',2)
hold on
plot(1:100,Mean_DP_Omni_Opt(:,2),'g:','LineWidth',2)
plot(1:100,Mean_DP_Omni_Opt(:,3),'r--','LineWidth',2)
plot(1:100,Mean_DP_Omni_Opt(:,4),'k-.','LineWidth',2);
hold off
axis([0 100 0.23 0.28]);
set(gca,'Fontname','Times New Roman','Fontsize',figure_FontSize)

P1=legend ('Region 1', 'Region 2' ,'Region 3' ,'Region 4');

 L_xoff   =  0.64    ;
 L_yoff   =  0.69    ;
 L_width  =  0.12  ;
 L_height =  0.001  ;
    
        set(P1,'Position',[L_xoff,L_yoff,L_width, L_height])
        set(P1,'FontName','Times New roman','FontSize',figure_FontSize-5 )


title ('Omni-opt','FontName','Times New roman','Fontsize',figure_FontSize)
xlabel ('Generation number','FontName','Times New roman','Fontsize',figure_FontSize)
ylabel('Proportion','FontName','Times New roman','Fontsize',figure_FontSize)

set(gcf,'Units','centimeters');%把图像大小单位用厘米衡量
set(gcf,'Position',[5 5  12.9 10.5]);
set(gca,'Position',[.15 .12 .81 .79]);

figure
plot(1:100,Mean_DP_Decision_niched_NSGAII(:,1),'b-','LineWidth',2)
hold on
plot(1:100,Mean_DP_Decision_niched_NSGAII(:,2),'g:','LineWidth',2)
plot(1:100,Mean_DP_Decision_niched_NSGAII(:,3),'r--','LineWidth',2)
plot(1:100,Mean_DP_Decision_niched_NSGAII(:,4),'k-.','LineWidth',2);
hold off
axis([0 100 0.23 0.28]);
set(gca,'Fontname','Times New Roman','Fontsize',figure_FontSize)
P1=legend ('Region 1', 'Region 2' ,'Region 3' ,'Region 4');

 L_xoff   =  0.64    ;
 L_yoff   =  0.69    ;
 L_width  =  0.15  ;
 L_height =  0.01  ;
    
        set(P1,'Position',[L_xoff,L_yoff,L_width, L_height])
        set(P1,'FontName','Times New roman','FontSize',figure_FontSize-5 )
title ('DN-NSGAII','FontName','Times New roman','Fontsize',figure_FontSize)
xlabel ('Generation number','FontName','Times New roman','Fontsize',figure_FontSize)
ylabel('Proportion','FontName','Times New roman','Fontsize',figure_FontSize)
set(gcf,'Units','centimeters');%把图像大小单位用厘米衡量
set(gcf,'Position',[5 5  12.9 10.5]);
set(gca,'Position',[.15 .12 .81 .79]);

% set(P1,'Box','Off' )
