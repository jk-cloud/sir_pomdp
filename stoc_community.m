function  [S,I,R,V,ISO,D,ICA,REWARD] = stoc_community(DNA)
%STOC_COMMUNITY stochastic simulation of a community
%   Driver routine for simulating the spread of a disease within a 
%   community of people using many realistations of the function
%
%    test_community
% 
%   to build averages.
%
%  usage:  S,I,R,V,ISO,D,ICA,REWARD] = test_stoc_community([DNA])
%
%      where DNA sets the reward matrix in the class Person
%      DNA con be used to interface the code using a genetic algorithm.
%      in order to optimze the reward matrix.s

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%
switch(nargin)
    case 1
      
    otherwise
      DNA=[];
end

M=12;      % realisations  
N=1000;   % PopulationSize
steps=30; % time steps

p=Person(0);
p.ConsistencyCheck;

%DNA=3*rand(1,p.GetNumberOfActions*p.GetNumberOfStates);
%DNA=round(DNA);

S=zeros(M,steps);   % susecptible
I=zeros(M,steps);   % infectious
R=zeros(M,steps);   % recovered
V=zeros(M,steps);   % vaccinated
ISO=zeros(M,steps); % isolated
D=zeros(M,steps);   % dead
ICA=zeros(M,steps); % intensive care
REWARD=zeros(M,steps); % total reward

if(license('test','Distrib_Computing_Toolbox'))
    parfor i=1:M
       [S(i,:),I(i,:),R(i,:),V(i,:),ISO(i,:),D(i,:),ICA(i,:),REWARD(i,:)]=community(N,steps,DNA);
    end
else
    for i=1:M
       [S(i,:),I(i,:),R(i,:),V(i,:),ISO(i,:),D(i,:),ICA(i,:),REWARD(i,:)]=community(N,steps,DNA);
    end
end


end





