function [ps,pf]=MO_PSO(func_name,VRmin,VRmax,Particle_Number,Max_Gen)
% MO_PSO: A multi-objective particle swarm optimization without ring topology and special crowding distance
%         each particle can exchange imformation with the whole population
%         in environment selection nondominated solutions are reserved
%       

%% Input£º
%                      Dimension                    Description
%      func_name       1 x length function name     the name of test function     
%      VRmin           1 x n_var                    low bound of decision variable
%      VRmax           1 x n_var                    up bound of decision variable
%      Particle_Number 1 x 1                        population size
%      Max_Gen         1 x 1                        maximum  generations

%% Output:
%                     Description
%      ps             Pareto set
%      pf             Pareto front
%%  Reference and Contact 
% Reference: [1]Caitong Yue, Boyang Qu and Jing Liang, "A Multi-objective Particle Swarm Optimizer Using Ring Topology for Solving Multimodal Multi-objective Problems",  IEEE Transactions on Evolutionary Computation, 2017.       
%            [2]Jing Liang, Caitong Yue, and Boyang Qu, ¡° Multimodal multi-objective optimization: A preliminary study¡±, IEEE Congress on Evolutionary Computation 2016, pp. 2454-2461, 2016.
% Contact: For any questions, please feel free to send email to zzuyuecaitong@163.com.
     

    Max_FES=Max_Gen*Particle_Number;
    n_PBA=5;              
    n_EXA=3*n_PBA;  %
    EXA=[];
    PBA{Particle_Number,1}=[];



    ps=Particle_Number;
    D=size(VRmin,2);
    cc=[2.05 2.05];   %acceleration constants
    iwt=0.7298;

    if length(VRmin)==1
        VRmin=repmat(VRmin,1,D);
        VRmax=repmat(VRmax,1,D);
    end
    mv=0.5*(VRmax-VRmin);
    VRmin=repmat(VRmin,ps,1);
    VRmax=repmat(VRmax,ps,1);
    Vmin=repmat(-mv,ps,1);
    Vmax=-Vmin;
    pos=VRmin+(VRmax-VRmin).*rand(ps,D);

    for i=1:ps;
    e(i,:)=feval(func_name,pos(i,:));
    end

    n_obj=size(e,2);
    particle=[pos,e];
    fitcount=ps;
    [D_index] = Non_dominated(particle(:,D+1:D+n_obj));
    EXA=particle(D_index,:);
for i=1:ps
    PBA{i,1}=particle(i,:);
end

vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles


for i=2:Max_gen

    for k=1:ps
   
   %% Randomly choose one particle in PBA{i,1} as pbest
    PBA_i=PBA{k,1};                  
    p_index=randperm(size(PBA_i,1),1); 
    pbest=PBA_i(p_index,:); 
   %% Randomly choose one particle in EXA as nbest
    n_index=randperm(size(EXA,1),1);
    nbest=EXA(n_index,:);
    
    aa(k,:)=cc(1).*rand(1,D).*(pbest(1,1:D)-pos(k,:))+cc(2).*rand(1,D).*(nbest(1,1:D)-pos(k,:));
    vel(k,:)=iwt.*vel(k,:)+aa(k,:); 
    vel(k,:)=(vel(k,:)>mv).*mv+(vel(k,:)<=mv).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv)).*(-mv)+(vel(k,:)>=(-mv)).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    pos(k,:)=((pos(k,:)>=VRmin(1,:))&(pos(k,:)<=VRmax(1,:))).*pos(k,:)...
        +(pos(k,:)<VRmin(1,:)).*(VRmin(1,:)+0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D))+(pos(k,:)>VRmax(1,:)).*(VRmax(1,:)-0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D));

    e(k,:)=feval(func_name,pos(k,1:D));
    fitcount=fitcount+1;
    particle(k,1:D+n_obj)=[pos(k,:),e(k,:)];
   %% Update PBA, if particle(k,:) is not dominated by any other solution in PBA{i,1}, put it in PBA{i,1}
   %   if the number of particles in PBA{i,1} is more than presetting number, delete the extra number of particles randomly
     PBA_i=[PBA_i;particle(k,:)];                     
     tempindex=Non_dominated(PBA_i(:,D+1:D+n_obj));
     tempPBAi=PBA_i(tempindex,:);
     size_tempPBAi=size(tempPBAi,1);                 
     if size_tempPBAi>n_PBA
         delet_PBA=size_tempPBAi-n_PBA;
         delet_PBAindex=randperm(size_tempPBAi,delet_PBA);
         tempPBAi(delet_PBAindex,:)=[];
     end
     PBA{k,1}=tempPBAi;
   
    end
    %% Update EXA
        pbest_A=cell2mat(PBA);
        tempEXA=pbest_A;                    
         tempindex=Non_dominated(tempEXA(:,D+1:D+n_obj));
         tempEXA=tempEXA(tempindex,:);
         size_tempEXA=size(tempEXA,1);                 
         if size_tempEXA>n_EXA
             delet_EXA=size_tempEXA-n_EXA;
             delet_EXAindex=randperm(size_tempEXA,delet_EXA);
             tempEXA(delet_EXAindex,:)=[];
         end
         EXA=tempEXA;
    
    

    if fitcount>=Max_FES
        break;
    end

end

   tempindex=Non_dominated(EXA(:,D+1:D+n_obj));
   EXA=EXA(tempindex,:);
   ps=EXA(:,1:D);
   pf=EXA(:,D+1:D+n_obj);


function [D_index] = Non_dominated(S)    

[m, n] = size(S);
D = [];
D_index=[];
for i=1:m
    [s, t] = size(D);
    flag1 = 0;
    for j=1:s
        flag2 = StrictDominate(D(j, :), S(i, :), n);
        if flag2~=0        
            flag1 = 1;     
        end
    end
    if flag1==0
        id = 1;
        for j=1:s
            flag2 = StrictDominate(S(i, :), D(id, :), n);
            if flag2==1
                D(id, :) = [];
                D_index(id, :)=[];
                id = id - 1;
            end
            id = id + 1;
        end
        D = [D; S(i, :)];
        D_index=[D_index;i];
    end
end

function [flag] = StrictDominate(a, b, n)

flag = 1;
count = 0;
for i=1:n
    if a(i)>b(i)
        flag = 0;
    end
    if a(i)==b(i)
        count = count + 1;
    end
end
if count==n
    flag = 2;
end
%Write by Caitong Yue 2016.3.16