function  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community(N,steps,DNA)
%TEST_COMMUNITY Simulating a community
%   Driver routine for simulating the spread of a disease within a 
%   community of people.
%   Each of community member is represented by an instance of the class
%   Person.
%   Most of the computation is done within "Person".
%   
%   This routine computes for each community member the probability of
%   getting infected by other menbers of the community using a simple 
%
%       p= 1-exp(-a* "number of contacts with infectious")
%
%   The parameter a can and should be adjusted. 
%   Another parameter for modelling the number ContactsPerDay which models
%   the contact network of a person
%
%   usage:  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community(N,steps,[DNA])

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%

switch(nargin)
    case 1
      steps=30;
      DNA=[];
      PLOT=false;
    case 2
      DNA=[];
      PLOT=false; 
    otherwise 
      DNA=[];
      PLOT=true;
end

N=1000;    % PopulationSize
steps=30;     % time steps

P0=Person(0);
P0.ConsistencyCheck;
na=P0.GetNumberOfActions;

DNA=3*rand(1,na*P0.GetNumberOfStates);
DNA=round(DNA);

S=zeros(1,steps);   % susecptible
I=zeros(1,steps);   % infectious
R=zeros(1,steps);   % recovered
V=zeros(1,steps);   % vaccinated
ISO=zeros(1,steps); % isolated
D=zeros(1,steps);   % dead
ICA=zeros(1,steps); % intensive care
REWARD=zeros(1,steps); % total reward
A=zeros(na,steps);  % actions
%initialise community


for i=1:N
    P{i}=Person(0.01,DNA);
end

% contact matrix
ContactsPerDay=10;
C=ceil(N*rand(N,ContactsPerDay));
p=0.1*ones(1,N);

% evolve
for t=1:steps
    % update

    for i=1:N
        P{i}.UpdatePerson(p(i));
    end

    %collect
    sus=0; % susceptible
    inf=0; % infectious
    rec=0; % recovered
    vac=0; % vaccinated
    iso=0; % isolated
    dea=0; % dead
    inc=0; % intensive care
    rew=0; % reward
    
    for i=1:N
        if(P{i}.IsSusceptible)
            sus=sus+1;
        end
        if(P{i}.IsInfectious)
            inf=inf+1;
        end
        if(P{i}.IsRecovered)
            rec=rec+1;
        end
        if(P{i}.IsVaccinated)
            vac=vac+1;
        end
        if(P{i}.IsInIsolation)
            iso=iso+1;
        end
        if(P{i}.IsDead)
            dea=dea+1;
        end
        if(P{i}.IsInIntensiveCare)
            inc=inc+1;
        end
        rew=rew+P{i}.reward;
        A(P{i}.a,t)=A(P{i}.a,t)+1;
    end
    
    S(t)=sus;   % susecptible
    I(t)=inf;   % infectious
    R(t)=rec;   % recovered
    V(t)=vac;   % vaccinated
    ISO(t)=iso; % isolated
    D(t)=dea;   % dead
    ICA(t)=inc; % intensive care
    REWARD(t)=rew; % total reward

    
    % update probabilities of getting infected
    for i=1:N
       p(i)=0;
       for j=1:ContactsPerDay
         if(P{C(i,j)}.SpreadsInfection)
            p(i)=p(i)+1;
         end
         p(i)=1-exp(-0.05*p(i)); % 30% probability of getting infected when meeting 7 infected per day
       end
    end
    
end

if(PLOT)
    phi=-20;
    psi=30;
    figure(2)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','actions evolution of a community')
    clf
    A(1,:)=0;
    bar3(A)
    xlabel('time [days]');
    ylabel('action');
    lbs=P0.GetActions;
    yticklabels(lbs);
    
    
    figure(3)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','evolution of a community')
    clf
    subplot(9,1,1)

    
    
    subplot(9,1,2)
    plot(S)
    
    ylabel('susceptible')
    
    
    subplot(9,1,3)
    plot(I)
    
    ylabel('infectious')
    
    subplot(9,1,4)
    plot(R)
    
    ylabel('recovered')
    
    subplot(9,1,5)
    plot(V)
    ylabel('vaccinated')
    
    subplot(9,1,6)
    plot(ISO)
    
    ylabel('isolated')
    
    subplot(9,1,7)
    plot(D)
    
    ylabel('dead')
    
    subplot(9,1,8)
    plot(ICA)
    
    ylabel('intensive care')
    
    subplot(9,1,9)
    plot(REWARD)
    xlabel('time [day]')
    ylabel('total reward')
    
end
end




