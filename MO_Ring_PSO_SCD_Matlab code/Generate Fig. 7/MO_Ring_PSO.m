function [ps,pf]=MO_Ring_PSO(func_name,VRmin,VRmax,Particle_Number,Max_Gen)
% MO_PSO_SCD: A multi-objective particle swarm optimization with only special crowding distance
%         each particle can exchange imformation with the whole population
%         in environment selection, the particles are sorted according to nondominated relationship and special crowding distance
%         
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
n_NBA=3*n_PBA;
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

row_of_cell=ones(1,ps);
col_of_cell=size(particle,2);
PBA=mat2cell(particle,row_of_cell,col_of_cell);
fitcount=ps;
for i=1:ps
     if i==1
         tempNBA=PBA{ps,:};
         tempNBA=[tempNBA;PBA{1,:}];
         tempNBA=[tempNBA;PBA{2,:}];
     elseif i==ps
         tempNBA=PBA{ps-1,:};
         tempNBA=[tempNBA;PBA{ps,:}];
         tempNBA=[tempNBA;PBA{1,:}];
     else
         tempNBA=PBA{i-1,:};
         tempNBA=[tempNBA;PBA{i,:}];
         tempNBA=[tempNBA;PBA{i+1,:}];
     end
     tempindex=Non_dominated(tempNBA(:,D+1:D+n_obj));
     NBA{i,1}=tempNBA(tempindex,:);
     clear  tempNBA tempindex
end

vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles


for i=2:Max_Gen
    

    for k=1:ps
   
    %% Randomly choose one particle in PBA{i,1} as pbest
    PBA_k=PBA{k,1};                
    p_index=randperm(size(PBA_k,1),1); 
    pbest=PBA_k(p_index,:); 
    %% Randomly choose one particle in EXA as nbest
    NBA_k=NBA{k,:};
    n_index=randperm(size(NBA_k,1),1);
    nbest=NBA_k(n_index,:);
    
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
     PBA_k=[PBA_k;particle(k,:)];                     
     tempindex=Non_dominated(PBA_k(:,D+1:D+n_obj));
     tempPBAi=PBA_k(tempindex,:);
     size_tempPBAi=size(tempPBAi,1);             
     if size_tempPBAi>n_PBA
         delet_PBA=size_tempPBAi-n_PBA;
         delet_PBAindex=randperm(size_tempPBAi,delet_PBA);
         tempPBAi(delet_PBAindex,:)=[];
     end
     PBA{k,1}=tempPBAi;
     
     %% Update EXA
     if k==1
         tempNBA=PBA{ps,:};
         tempNBA=[tempNBA;PBA{1,:}];
         tempNBA=[tempNBA;PBA{2,:}];
     elseif k==ps
         tempNBA=PBA{ps-1,:};
         tempNBA=[tempNBA;PBA{ps,:}];
         tempNBA=[tempNBA;PBA{1,:}];
     else
         tempNBA=PBA{k-1,:};
         tempNBA=[tempNBA;PBA{k,:}];
         tempNBA=[tempNBA;PBA{k+1,:}];
     end
     tempindex=Non_dominated(tempNBA(:,D+1:D+n_obj));
     if length(tempindex)>=n_NBA
         NBA{k,1}=tempNBA(tempindex(1:n_NBA),:);
     else
          NBA{k,1}=tempNBA(tempindex,:);
     end
         
    
     
    end

 
 
    if fitcount>=Max_FES
        break;
    end

end
   EXA=cell2mat(NBA);
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