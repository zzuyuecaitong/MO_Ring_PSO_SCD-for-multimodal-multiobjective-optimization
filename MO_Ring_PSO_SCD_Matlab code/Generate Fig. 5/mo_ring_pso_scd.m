function [ps,pf,Dis_proportion]=mo_ring_pso_scd(func_name,VRmin,VRmax,Particle_Number,Max_Gen)
% MO_Ring_PSO_SCD: A multi-objective particle swarm optimization using ring topology for solving multimodal multi-objective optimization problems 
% Dimension: n_var --- dimensions of decision space
%            n_obj --- dimensions of objective space
%% Input:
%                      Dimension                    Description
%      func_name       1 x length function name     the name of test function     
%      VRmin           1 x n_var                    low bound of decision variable
%      VRmax           1 x n_var                    up bound of decision variable
%      n_obj           1 x 1                        dimensions of objective space
%      Particle_Number 1 x 1                        population size
%      Max_Gen         1 x 1                        maximum  generations

%% Output:
%                      Dimension                              Description
%      ps              population_size x n_var                Pareto set
%      pf              population_size x n_obj                Pareto front
%      Dis_proportion  num_generation  x num_of_regions       the proportion of population in different regions
%%  Reference and Contact 
% Reference: [1]Caitong Yue, Boyang Qu and Jing Liang, "A Multi-objective Particle Swarm Optimizer Using Ring Topology for Solving Multimodal Multi-objective Problems",  IEEE Transactions on Evolutionary Computation, 2017.       
%            [2]Jing Liang, Caitong Yue, and Boyang Qu, ¡° Multimodal multi-objective optimization: A preliminary study¡±, IEEE Congress on Evolutionary Computation 2016, pp. 2454-2461, 2016.
% Contact: For any questions, please feel free to send email to zzuyuecaitong@163.com.       

D=size(VRmin,2);
Max_FES=Max_Gen*Particle_Number;
n_PBA=5;                  
n_NBA=3*n_PBA;
PBA{Particle_Number,1}=[];


cc=[2.05 2.05];   %acceleration constants
iwt=0.7298;

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,Particle_Number,1);
VRmax=repmat(VRmax,Particle_Number,1);
Vmin=repmat(-mv,Particle_Number,1);
Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(Particle_Number,D);

for i=1:Particle_Number;
e(i,:)=feval(func_name,pos(i,:));
end

n_obj=size(e,2);
particle=[pos,e];

row_of_cell=ones(1,Particle_Number);
col_of_cell=size(particle,2);
PBA=mat2cell(particle,row_of_cell,col_of_cell);
NBA=PBA;
fitcount=Particle_Number;
vel=Vmin+2.*Vmax.*rand(Particle_Number,D);%initialize the velocity of the particles
EXA=[];
tempEXA=[];

for i=1:Max_Gen
    
    %% Update NBA
    for j=1:Particle_Number
     if j==1
         tempNBA=PBA{Particle_Number,:};
         tempNBA=[tempNBA;PBA{1,:}];
         tempNBA=[tempNBA;PBA{2,:}];
     elseif j==Particle_Number
         tempNBA=PBA{Particle_Number-1,:};
         tempNBA=[tempNBA;PBA{Particle_Number,:}];
         tempNBA=[tempNBA;PBA{1,:}];
     else
         tempNBA=PBA{j-1,:};
         tempNBA=[tempNBA;PBA{j,:}];
         tempNBA=[tempNBA;PBA{j+1,:}];
     end
     NBA_j=NBA{j,1};
     tempNBA=[tempNBA;NBA_j];
     tempNBA=non_domination_sort_crowd_dist(tempNBA(:,1:D+n_obj), n_obj, D);
     if size(tempNBA,1)>n_NBA
         NBA{j,1}=tempNBA(1:n_NBA,:);
     else
        NBA{j,1}=tempNBA;
     end
          
    end
    
%  disp(['mo_ring_pso_Omni_4y21r  gen=' num2str(i)  ])
Pbest=[];
Nbest=[];
    for k=1:Particle_Number
   
    
    PBA_k=PBA{k,1};                  
    pbest=PBA_k(1,:); 
    Pbest=[Pbest;pbest];

    NBA_k=NBA{k,:};
    nbest=NBA_k(1,:);
    Nbest=[Nbest;nbest];
    
    Population=cell2mat(NBA);                    
   
    
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
     %% Update PBA
     PBA_k=[PBA_k(:,1:D+n_obj);particle(k,:)];                    
     PBA_k = non_domination_sort_crowd_dist(PBA_k(:,1:D+n_obj), n_obj, D);
     if size(PBA_k,1)>n_PBA
         PBA{k,1}=PBA_k(1:n_PBA,:);
     else
         PBA{k,1}=PBA_k;
     end
     

    end
%     figure(1)
%     plot(Population(:,1),Population(:,2),'bo',Pbest(:,1),Pbest(:,2),'g+',Nbest(:,1),Nbest(:,2),'r*');
%     legend ('population','pbest','nbest')
%     title ([ num2str(i) '^th Generation'])
 Dis_proportion(i,:)=calculate_region_proportion (Population(:,1:D));% num_generation x num_of_regions
  
    if fitcount>Max_FES
        break;
    end

end

    tempEXA=cell2mat(NBA);                    
   
    tempEXA=[tempEXA;EXA];                    
    tempEXA=non_domination_sort_crowd_dist(tempEXA(:,1:D+n_obj), n_obj, D);
     if size(tempEXA,1)>Particle_Number
         EXA=tempEXA(1:Particle_Number,:);
     else
        EXA=tempEXA;
     end

ps=EXA(:,1:D);
pf=EXA(:,D+1:D+n_obj);
end

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
end
function f = non_domination_sort_crowd_dist(x, M, V)


%%
[N, m] = size(x);
clear m

% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = [];

%% Non-Dominated sort. 

for i = 1 : N
    % Number of individuals that dominate this individual
    individual(i).n = 0; 
    % Individuals which this individual dominate
    individual(i).p = [];
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k))
                dom_less = dom_less + 1;
            elseif (x(i,V + k) == x(j,V + k))  
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).p = [individual(i).p j];
        end
    end   
    if individual(i).n == 0
        x(i,M + V + 1) = 1;
        F(front).f = [F(front).f i];
    end
end
% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)
       if ~isempty(individual(F(front).f(i)).p)
        	for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;

%% Crowding distance

for front = 1 : (length(F) - 1)
%    objective = [];
    crowd_dist_obj = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;
    % Sort each individual based on the objective
    sorted_based_on_objective = [];
    for i = 1 : M+V
        [sorted_based_on_objective, index_of_objectives] = ...
            sort(y(:,i));
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
        end
        f_max = ...
            sorted_based_on_objective(length(index_of_objectives), i);
        f_min = sorted_based_on_objective(1,  i);

        if length(index_of_objectives)==1
            y(index_of_objectives(1),M + V + 1 + i) = 1;
          
        elseif i>V

            y(index_of_objectives(1),M + V + 1 + i) = 1;
            y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)=0;
        else
       
             y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...
                = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
            sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min);
             y(index_of_objectives(1),M + V + 1 + i)=2*(sorted_based_on_objective(2, i)-...
            sorted_based_on_objective(1, i))/(f_max - f_min);
        end
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1, i);
            previous_obj  = sorted_based_on_objective(j - 1,i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + V + 1 + i) = 1;
            else
                y(index_of_objectives(j),M + V + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    %Calculate distance in x space
    crowd_dist_var = [];
    crowd_dist_var(:,1) = zeros(length(F(front).f),1);
    for i = 1 : V
        crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,M + V + 1 + i);
    end
    crowd_dist_var=crowd_dist_var./V;
    avg_crowd_dist_var=mean(crowd_dist_var);
    %Calculate distance in f space
    crowd_dist_obj = [];
    crowd_dist_obj(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,M + V + 1+V + i);
    end
    crowd_dist_obj=crowd_dist_obj./M;
    avg_crowd_dist_obj=mean(crowd_dist_obj);
    crowd_dist=zeros(length(F(front).f),1);
    for i = 1 : length(F(front).f)
        if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
            crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i));
        else
            crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i));
        end
    end
    y(:,M + V + 2) = crowd_dist;
    y(:,M+V+3)=crowd_dist_var;
    y(:,M+V+4)=crowd_dist_obj;
    [~,index_sorted_based_crowddist]=sort(crowd_dist,'descend');
    y=y(index_sorted_based_crowddist,:);
    y = y(:,1 : M + V+4 );
    z(previous_index:current_index,:) = y;
end
f = z();
end

%Write by Caitong Yue 2016.3.16