classdef community  < handle
%  [S,I,R,V,ISO,D,ICA,REWARD,A] = community(N,steps,DNA)
%COMMUNITY Simulating a community
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


properties
    Name = []
end


properties(SetAccess=protected)
    N=[]; % PopulationSize
    M=[]; % number of simulation realizations
    steps=[]; % time steps (days)
    DNA=[];
    S=[];   % susecptible
    I=[];   % infectious
    R=[];   % recovered
    V=[];   % vaccinated
    ISO=[]; % isolated
    D=[];   % dead
    ICA=[]; % intensive care
    REWARD=[]; % total reward
    A=[];  % actions
    
    P=[]; % list of persons
    
    C=[]; % contact matrix
    ContactsPerDay =[];
    ages =[];  % ages of persons
    NeedForIntensiveCareAge=[];     % interpolation vectors for intensive care
    NeedForIntensiveCare=[];    
    IntensiveCareTimeAge=[];
    IntensiveCareTimeDays=[];
    IntensiveCareRecoveryAge=[];
    IntensiveCareRecoveryProb=[];
end



methods
    %%
    function obj = community(DNA)
        %Community Construct an instance of this class
        %s
        switch nargin
            case 1
                obj.DNA=DNA;

            otherwise

        end
        obj.Name = "Community";
        obj.Properties;
        obj.M=1;
    end
    %%
    function obj = Properties(obj)
        obj.N=1000;
        obj.steps=30;
        obj.ContactsPerDay = 10;
    end
    %%
    
    function obj = SetPopulationSize(obj,N)
        obj.N = N;
    end
    
    %
    function obj = SetSimulationSteps(obj,steps)
        obj.steps = steps;
    end  
    %
    function obj = SetNumberOfRealizations(obj,M)
        obj.M = M;
    end  
    %
    function obj = Initialize(obj)
        P0=Person(0);
        P0.ConsistencyCheck;
        na=P0.GetNumberOfActions;

        obj.S=zeros(obj.M,obj.steps);   % susecptible
        obj.I=zeros(obj.M,obj.steps);   % infectious
        obj.R=zeros(obj.M,obj.steps);   % recovered
        obj.V=zeros(obj.M,obj.steps);   % vaccinated
        obj.ISO=zeros(obj.M,obj.steps); % isolated
        obj.D=zeros(obj.M,obj.steps);   % dead
        obj.ICA=zeros(obj.M,obj.steps); % intensive care
        obj.REWARD=zeros(obj.M,obj.steps); % total reward
        obj.A=zeros(na,obj.steps);  % actions
        
        
        
        %initialise community
        obj.InitializeAges;
        for i=1:obj.N
            obj.P{i}=Person(0.01,obj.DNA,obj.ages(i));
            obj.P{i}.SetNeedForIntensiveCareBasedOnAge(obj.NeedForIntensiveCareAge,obj.NeedForIntensiveCare);
            obj.P{i}.SetIntensiveCareRecoveryBasedOnAge(obj.IntensiveCareRecoveryAge,obj.IntensiveCareRecoveryProb);
            obj.P{i}.SetIntensiveCareTimeBasedOnAge(obj.IntensiveCareTimeAge,obj.IntensiveCareTimeDays);
            obj.P{i}.SetParametersBasedOnAge;  
        end
        
        
        % contact matrix
        obj.ContactsPerDay=10;
        obj.C=ceil(obj.N*rand(obj.N,obj.ContactsPerDay));
        
        
    end
    
    %
    function obj = InitializeAges(obj)
       X = randn(1,obj.N);
       X = X + min(X);
       
       obj.ages =  50*ones(1,obj.N);
    end
    %
    function obj = PlotAges(obj)
       plot(obj.ages) 
    end
    %
    function obj = Evolve(obj)
        P0=Person(0);
        na=P0.GetNumberOfActions;
        % Evolve uses Evolution and assigns enough memory
        % for the results.
        % The results are copied back to class properties.
        % This is clumsy but nececcessary due to the usage of the 
        % parfor loop.
        M=obj.M;
        S=zeros(M,obj.steps);   % susecptible
        I=zeros(M,obj.steps);   % infectious
        R=zeros(M,obj.steps);   % recovered
        V=zeros(M,obj.steps);   % vaccinated
        ISO=zeros(M,obj.steps); % isolated
        D=zeros(M,obj.steps);   % dead
        ICA=zeros(M,obj.steps); % intensive care
        REWARD=zeros(M,obj.steps); % total reward
        A=zeros(na,obj.steps);  % actions

        if(M>1 && license('test','Distrib_Computing_Toolbox'))
            parfor i=1:M
                [S(i,:),I(i,:),R(i,:),V(i,:),ISO(i,:),D(i,:),ICA(i,:),REWARD(i,:),A]=obj.Evolution;
            end
        else
            for i=1:M
                [S(i,:),I(i,:),R(i,:),V(i,:),ISO(i,:),D(i,:),ICA(i,:),REWARD(i,:),A]=obj.Evolution;
            end
        end
        
        obj.S = S;
        obj.I = I;
        obj.R = R;
        obj.V =V;
        obj.ISO = ISO;
        obj.D =D;
        objICA = ICA;
        obj.REWARD =REWARD;
        obj.A =A;
    end
    
    function [S,I,R,V,ISO,D,ICA,REWARD,A] = ReturnResults(obj)
       
        S = obj.S;
        I = obj.I;
        R = obj.R;
        V = obj.V;
        ISO= obj.ISO;
        D = obj.D;
        ICA =obj.ICA;
        REWARD = obj.REWARD;
        A = obj.A;
    end
        
    %
    function [S,I,R,V,ISO,D,ICA,REWARD,A] = Evolution(obj)
        P0=Person(0);
        na=P0.GetNumberOfActions;
        
        S=zeros(1,obj.steps);   % susecptible
        I=zeros(1,obj.steps);   % infectious
        R=zeros(1,obj.steps);   % recovered
        V=zeros(1,obj.steps);   % vaccinated
        ISO=zeros(1,obj.steps); % isolated
        D=zeros(1,obj.steps);   % dead
        ICA=zeros(1,obj.steps); % intensive care
        REWARD=zeros(1,obj.steps); % total reward
        A=zeros(na,obj.steps);  % actions
        p=0.1*ones(1,obj.N);
        for t=1:obj.steps
            % update
            
            for i=1:obj.N
                obj.P{i}.UpdatePerson(p(i));
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
            
            for i=1:obj.N
                if(obj.P{i}.IsSusceptible)
                    sus=sus+1;
                end
                if(obj.P{i}.IsInfectious)
                    inf=inf+1;
                end
                if(obj.P{i}.IsRecovered)
                    rec=rec+1;
                end
                if(obj.P{i}.IsVaccinated)
                    vac=vac+1;
                end
                if(obj.P{i}.IsInIsolation)
                    iso=iso+1;
                end
                if(obj.P{i}.IsDead)
                    dea=dea+1;
                end
                if(obj.P{i}.IsInIntensiveCare)
                    inc=inc+1;
                end
                rew=rew+obj.P{i}.reward;
                obj.A(obj.P{i}.a,t)=obj.A(obj.P{i}.a,t)+1;
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
            for i=1:obj.N
                p(i)=0;
                for j=1:obj.ContactsPerDay
                    if(obj.P{obj.C(i,j)}.SpreadsInfection)
                        p(i)=p(i)+1;
                    end
                    p(i)=1-exp(-0.5*p(i)); % 30% probability of getting infected when meeting 7 infected per day
                end
            end
            
        end
    end
    %%
    function obj=SetNeedForIntensiveCareBasedOnAge(obj,x,y)
        obj.NeedForIntensiveCareAge = x;
        obj.NeedForIntensiveCare    = y;
    end
    %%
    function obj=SetIntensiveCareRecoveryBasedOnAge(obj,x,y)
        obj.IntensiveCareRecoveryAge  = x;
        obj.IntensiveCareRecoveryProb = y;
    end
    %%
    function obj=SetIntensiveCareTimeBasedOnAge(obj,x,y)
        obj.IntensiveCareTimeAge  = x;
        obj.IntensiveCareTimeDays = y;
    end
    %%
    
end

end




