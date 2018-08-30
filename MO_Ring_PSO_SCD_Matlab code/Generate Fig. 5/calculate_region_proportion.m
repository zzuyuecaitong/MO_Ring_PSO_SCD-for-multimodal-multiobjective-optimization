function [Proportion]=calculate_region_proportion (Population)

%calculate_region_proportion: Calculate the proportion of population in different regions
%% Input:
%                      Dimension                    Description
%      Population      population_size x n_var      the name of test function     


%% Output:
%                   Dimension                    Description  
%      Proportion   1 x num_of_regions           the proportion of solutions in each region
  
bound1=[-1 0
        1  2];%bound1: n_var x 2;  the minimum is put in the first column, the maximum is put in the second column,region1 -1<x1<0,1<x2<2
bound2=[0 1
        1 2];
bound3=[-1 0
        0 1];
bound4=[0 1
        0 1];
 n_region=4;
count_r=zeros(1,n_region);%count_r: 1 x n_region£¬the number of solutions in the first region is put in the first column of count_r, and so on.
for i=1:size(Population,1)
    if Population(i,1)>=bound1(1,1)&&Population(i,1)<=bound1(1,2)&&Population(i,2)>=bound1(2,1)&&Population(i,2)<=bound1(2,2)
        count_r(1,1)=count_r(1,1)+1;%the first region count + 1
    elseif Population(i,1)>bound2(1,1)&&Population(i,1)<=bound2(1,2)&&Population(i,2)>=bound2(2,1)&&Population(i,2)<=bound2(2,2)
        count_r(1,2)=count_r(1,2)+1;%the second region count + 1
    elseif Population(i,1)>=bound3(1,1)&&Population(i,1)<=bound3(1,2)&&Population(i,2)>=bound3(2,1)&&Population(i,2)<bound3(2,2)
        count_r(1,3)=count_r(1,3)+1;
    elseif Population(i,1)>bound4(1,1)&&Population(i,1)<=bound4(1,2)&&Population(i,2)>=bound4(2,1)&&Population(i,2)<bound4(2,2)
        count_r(1,4)=count_r(1,4)+1;
    end
end
Proportion=count_r./sum(count_r);
end
