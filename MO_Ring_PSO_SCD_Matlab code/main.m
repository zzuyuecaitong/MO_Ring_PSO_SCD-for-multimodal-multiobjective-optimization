%% Add path
addpath(genpath('MM_testfunctions/'));
addpath(genpath('Indicator_calculation/'));
 clear all
  clc
%   rand('state',sum(100*clock));
  global fname
  N_function=11;% number of test function
  popsize=800;
  Max_evaluation=80000;
  Max_Gen=fix(Max_evaluation/popsize);
  % Note: It may take a long time to run all 11 test functions and with
  % population size 800 and generation 100. You can change N_function to 1,
  %  popsize to 100, Max_evaluation to 1000, to see how the MO_Ring_PSO_SCD
  %  works.
 
 
 for i=1:N_function
    switch i
        case 1
            fname='MMF1';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
            repoint=[2,2]; % reference point used to calculate the Hypervolume
            load('MMF1truePSPF.mat');
        case 2
            fname='MMF2';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 2];
            repoint=[2,2];
            load('MMF2truePSPF.mat');
        case 3
            fname='MMF3';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 1.5];
            repoint=[2,2];
            load('MMF3truePSPF.mat');
        case 4
            fname='MMF4';
            n_obj=2;
            n_var=2;
            xl=[-1 0];
            xu=[1 2];
            repoint=[2,2];
            load('MMF4truePSPF.mat');
        case 5
            fname='MMF5';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 3];
            repoint=[2,2];
            load('MMF5truePSPF.mat');
         case 6
            fname='MMF6';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 2];
            repoint=[2,2];
            load('MMF6truePSPF.mat');
        case 7
            fname='MMF7';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 1];
            repoint=[2,2];
            load('MMF7truePSPF.mat');
         case 8
            fname='MMF8';
            n_obj=2;
            n_var=2;
            xl=[-pi 0];
            xu=[pi 9];
            repoint=[2,2];
            load('MMF8truePSPF.mat');
         case 9
            fname='SYM_PART_simple';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[2,2];
            load('SYM_PART_simple_turePSPF.mat');
         case 10
            fname='SYM_PART_rotated';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[2,2];
            load('SYM_PART_rotatedtruePSPF.mat');
        case 11
            fname='Omni_test';
            n_obj=2;
            n_var=3;
            xl=[0 0 0];
            xu=[6 6 6];
            repoint=[5,5];
            load('Omni_testtruePSPF.mat');
    end
   fprintf('Running test function: %s \n', fname);
   %% Search the PSs using MO_Ring_PSO_SCD
    [ps,pf]=MO_Ring_PSO_SCD(fname,xl,xu,n_obj,popsize,Max_Gen);
   %% Indicators
     hyp=Hypervolume_calculation(pf,repoint);
     IGDx=IGD_calculation(ps,PS);
     CR=CR_calculation(ps,PS);
     PSP=CR/IGDx;% Eq. (8) in the paper
   
   %% Plot figure
       figure
       plot(ps(:,1),ps(:,2),'o');
       hold on;
       plot(PS(:,1),PS(:,2),'r+');
       legend 'Obtained PS' 'True PS'
       title (fname);
end


