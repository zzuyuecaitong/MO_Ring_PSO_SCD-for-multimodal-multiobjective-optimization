
 clear all
  clc
%   rand('state',sum(100*clock));
  global fname

  popsize=800;
  Max_Gen=fix(80000/popsize);
 
  %% Take MMF3 for example
            fname='MMF3';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 1.5];
            repoint=[2,2];
            load('MMF3truePSPF.mat');
       
   fprintf('Running test function: %s \n It will take period of time... \n', fname);
   %% Search the PSs using MO_Ring_PSO_SCD
    [ps,pf]=MO_Ring_PSO_SCD(fname,xl,xu,n_obj,popsize,Max_Gen);
    WHPS.MO_Ring_PSO_SCD.MMF3=ps;
    WHPF.MO_Ring_PSO_SCD.MMF3=pf;
   %% Search the PSs using MO_PSO_SCD
    [ps,pf]=MO_PSO_SCD(fname,xl,xu,popsize,Max_Gen);
    WHPS.MO_PSO_SCD.MMF3=ps;
    WHPF.MO_PSO_SCD.MMF3=pf;
   %% Search the PSs using MO_Ring_PSO
    [ps,pf]=MO_Ring_PSO(fname,xl,xu,popsize,Max_Gen);
    WHPS.MO_Ring_PSO.MMF3=ps;
    WHPF.MO_Ring_PSO.MMF3=pf;
   %% Search the PSs using MO_PSO 
    [ps,pf]=MO_PSO(fname,xl,xu,popsize,Max_Gen);
    WHPS.MO_PSO.MMF3=ps;
    WHPF.MO_PSO.MMF3=pf;
  
 save Fig7data WHPS WHPF  
 %% Plot Fig. 7 (a)-(d)
  figure(1)
       plot(WHPS.MO_Ring_PSO_SCD.MMF3(:,1),WHPS.MO_Ring_PSO_SCD.MMF3(:,2),'*');
       xlabel ('{\itx}_1','FontName','Times New Roman','FontSize',10);  
       ylabel ('{\itx}_2','FontName','Times New Roman','FontSize',10);  
       set(get(gca,'YLabel'),'Rotation', -pi/2);
  figure(2)
       plot(WHPS.MO_PSO_SCD.MMF3(:,1),WHPS.MO_PSO_SCD.MMF3(:,2),'*');
       xlabel ('{\itx}_1','FontName','Times New Roman','FontSize',10);  
       ylabel ('{\itx}_2','FontName','Times New Roman','FontSize',10);  
       set(get(gca,'YLabel'),'Rotation', -pi/2);
  figure(3)
       plot(WHPS.MO_Ring_PSO.MMF3(:,1),WHPS.MO_Ring_PSO.MMF3(:,2),'*');
       xlabel ('{\itx}_1','FontName','Times New Roman','FontSize',10);  
       ylabel ('{\itx}_2','FontName','Times New Roman','FontSize',10);  
       set(get(gca,'YLabel'),'Rotation', -pi/2);
 figure(4)
       plot(WHPS.MO_PSO.MMF3(:,1),WHPS.MO_PSO.MMF3(:,2),'*');
       xlabel ('{\itx}_1','FontName','Times New Roman','FontSize',10);  
       ylabel ('{\itx}_2','FontName','Times New Roman','FontSize',10);  
       set(get(gca,'YLabel'),'Rotation', -pi/2);

