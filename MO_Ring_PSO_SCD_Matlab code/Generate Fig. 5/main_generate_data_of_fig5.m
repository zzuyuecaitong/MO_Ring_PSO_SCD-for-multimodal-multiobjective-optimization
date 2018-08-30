  clear all
  clc
  rand('state',sum(100*clock));
  global fname
  N_function=7;
  runtimes=20;
  popsize=800;
  Max_Gen=fix(80000/popsize);

        fname='MMF4';
        n_obj=2;
        n=2;
        xl=[-1 0];
        xu=[1 2];
        repoint=[2,2];


    for j=1:runtimes
     %%  
     fprintf(2,' runtimes= %u/%u\n',j,runtimes)

    
      %%  mo_ring_pso_scd
     [ps,pf,DP_mo_ring_pso_scd(:,:,j)]=mo_ring_pso_scd(fname,xl,xu,popsize,Max_Gen);%DP_mo_ring_pso_scd行是代数，列是区域编号
     WHPS.Mo_ring_pso_scd.testf5=ps;
     
      %% Omni_Opt
     [ps,pf,DP_Omni_Opt(:,:,j)]=Omni_Opt(fname,n_obj,n,xl,xu,popsize,Max_Gen);
     WHPS.Omni_Opt.testf5=ps;
     
     %% Decision_niched_NSGAII
     [ps,pf,DP_Decision_niched_NSGAII(:,:,j)]=Decision_niched_NSGAII(fname,n_obj,n,xl,xu,popsize,Max_Gen);
    
     WHPS.Decision_niched_NSGAII.testf5=ps;
       
    end


save result_distribution_of_threealgorithms_20times DP_mo_ring_pso_scd DP_Omni_Opt DP_Decision_niched_NSGAII 
% generate figure
Mean_DP_mo_ring_pso_scd=mean(DP_mo_ring_pso_scd,3);
Mean_DP_Omni_Opt=mean(DP_Omni_Opt,3);
Mean_DP_Decision_niched_NSGAII=mean(DP_Decision_niched_NSGAII,3);
figure
plot(1:100,Mean_DP_mo_ring_pso_scd(:,1),1:100,Mean_DP_mo_ring_pso_scd(:,2),':',1:100,Mean_DP_mo_ring_pso_scd(:,3),'--',1:100,Mean_DP_mo_ring_pso_scd(:,4),'k-.');
axis([0 100 0.23 0.28]);
P1=legend ('Region1', 'Region2' ,'Region3' ,'Region4');
set(P1,'FontName','Times New roman')
title ('MO-Ring-PSO-SCD','FontName','Times New roman')
xlabel ('Generation','FontName','Times New roman')
ylabel('Proportion','FontName','Times New roman')

figure
plot(1:100,Mean_DP_Omni_Opt(:,1),1:100,Mean_DP_Omni_Opt(:,2),':',1:100,Mean_DP_Omni_Opt(:,3),'--',1:100,Mean_DP_Omni_Opt(:,4),'k-.');
axis([0 100 0.23 0.28]);
P1=legend ('Region1', 'Region2' ,'Region3' ,'Region4');
set(P1,'FontName','Times New roman')
title ('Omni-opt','FontName','Times New roman')
xlabel ('Generation','FontName','Times New roman')
ylabel('Proportion','FontName','Times New roman')

figure
plot(1:100,Mean_DP_Decision_niched_NSGAII(:,1),1:100,Mean_DP_Decision_niched_NSGAII(:,2),':',1:100,Mean_DP_Decision_niched_NSGAII(:,3),'--',1:100,Mean_DP_Decision_niched_NSGAII(:,4),'k-.');
axis([0 100 0.23 0.28]);
P1=legend ('Region1', 'Region2' ,'Region3' ,'Region4');
set(P1,'FontName','Times New roman')
title ('DN-NSGAII','FontName','Times New roman')
xlabel ('Generation','FontName','Times New roman')
ylabel('Proportion','FontName','Times New roman')




