classdef Person  < handle
    %PERSON class for POMDP
    %   The Person Class can be used simulate the state transition
    %   of a single person
    %   Initialization
    %
    %   P= Person(prob of getting infected)
    %   P= Person(0.1)
    %
    %
    %   available states
    %   P.DisplayTransitionIndizes
    %
    %   available observations
    %   P.DisplayObservationIndizes
    %
    %   available actions
    %   P.DisplayActions
    %
    
    %
    %  (c) 2020 Jens Kappey and the sir_pomdp contributors.
    %
    properties
        Name = []
    end
    %%
    properties(SetAccess=protected)
        notdetect=.999;                  % probability of not detecting the infectious state correctly without a test set [0.0 1.0]
        notdetectTestSet=0.2;           % probability of not detecting the infectious state correctly with a test set [0.0 1.0]
        IncubationPeriod=5;             % number of days of after infection, when symptoms show   Covid-19:  mean 5 days 2 to 20 days
        DiseaseDuration=10;             % duration in days of having no pathogens after infection
        vacc=0.0;                       % probability of successfully vaccinating a person [0.0 1.0]
        IntensiveCare=0.01;             % probability of needing intensive care when infected [0.0 1.0]
        IntensiveCareTime=10;           % time in days needing intensive care after being taken in intensive care
        IntensiveCareRecovery=0.995;     % probability of recovering when in intensive care  [0.0 1.0]
        IsolationDuration=10;           % time in days for someone being isolated at home ill or not
        notisolationSus=0.0;               % probability of not isolating oneself when susceptible
        notisolationInf=0.0;               % probability of not isolating oneself when infectious ... TODO not implemented yet
        
        ns=[];
        nib=[];
        nie=[];
        nr=[];
        nv=[];
        nitb=[];
        nite=[];
        nisob=[];
        nisoe=[];
        nisosusb=[];
        nisosuse=[];
        nd=[];
        n=[];
        
        ons=[];  % susceptible
        oni=[];  % infected
        onic=[]; % intensive care
        onv=[];  % vaccinated
        ond=[];  % dead
        onr=[];  % recovered
        oniso=[];% isolated
        no=[];   % total
        
        na=[];   % number of possible actions
        
        s=[];    % state
        z=[];    % observation
        b=[];    % belief
        a=[];    % action
        
        p=[];    % probability of getting infected by community
        p0=[];   % probability of getting infected from outside of community
        policy=[]; % policy of taken an action given state, observation, belief about state
        
        T=[];   % transition probabilities
        TR=[];  % transition realisation
        
        O=[];    % observation probabilities
        OR=[];   % observation realisation
        
        R=[];    % reward matrix
        reward=[]; % cumulative reward
        
        DNA=[];
        
        vectorize=false; % should be set to false. The alternative code is
        % much slower
    end
    
    methods
        %%
        function obj = Person(p,DNA)
            %PERSON Construct an instance of this class
            %
            switch nargin
                case 1
                    obj.p = p;
                    
                case 2
                    obj.p = p;
                    obj.DNA=DNA;
                    
                otherwise
                    
                    obj.p=0.01;
                    
            end
            obj.Name = "Person";
            obj.Properties;
        end
        %%
        function obj = Properties(obj)
            obj.ns=1;
            obj.nib=2;
            obj.nie=obj.nib+obj.DiseaseDuration;
            obj.nr=obj.nie+1;
            obj.nv=obj.nr+1;
            obj.nitb=obj.nv+1;
            obj.nite=obj.nitb+obj.IntensiveCareTime;
            obj.nisob=obj.nite+1;
            obj.nisoe=obj.nisob+obj.IsolationDuration;
            obj.nisosusb=obj.nisoe+1;
            obj.nisosuse=obj.nisosusb+obj.IsolationDuration;
            obj.nd=obj.nisosuse+1;
            obj.n=obj.nd;
            
            obj.ons=1;      % susceptible
            obj.oni=obj.ons+1;  % infected
            obj.onic=obj.oni+1; % intensive care
            obj.onv=obj.onic+1; % vaccinated
            obj.ond=obj.onv+1;  % dead
            obj.onr=obj.ond+1;  % recovered
            obj.oniso=obj.onr+1;  % isolated
            obj.no=obj.oniso;
            
            obj.s=zeros(obj.n,1);
            obj.s(obj.ns)=1;
            
            obj.na=5;
            
            obj.z=zeros(obj.no,1);
            obj.z(obj.ons)=1;
            obj.b=zeros(size(obj.s)); obj.b(1)=1; %ones(size(obj.s))/length(obj.s);
            obj.a=1; % do nothing as initial action;
            obj.reward=0;
            obj.InitRewardMatrix;
        end
        %%
        function str = DisplayTransitionIndizes(obj)
            i=1;
            str{i}=sprintf('susceptible: %i',obj.ns);i=i+1;
            str{i}=sprintf('infectious begin: %i',obj.nib);i=i+1;
            str{i}=sprintf('infectious end: %i',obj.nie);i=i+1;
            str{i}=sprintf('recovered: %i',obj.nr);i=i+1;
            str{i}=sprintf('vaccinated: %i',obj.nv);i=i+1;
            str{i}=sprintf('intensive care begin: %i',obj.nitb);i=i+1;
            str{i}=sprintf('intensive care end: %i',obj.nite);i=i+1;
            str{i}=sprintf('isolation infectious begin: %i',obj.nisob);i=i+1;
            str{i}=sprintf('isolation infectious end: %i',obj.nisoe);i=i+1;
            str{i}=sprintf('isolation susceptible begin: %i',obj.nisosusb);i=i+1;
            str{i}=sprintf('isolation susceptible end: %i',obj.nisosuse);i=i+1;
            str{i}=sprintf('dead: %i',obj.nd);i=i+1;
            for k=1:i-1
                disp(str{k})
            end
            
        end
        %%
        function str = DisplayObservationIndizes(obj)
            i=1;
            str{i}=sprintf('obs susceptible: %i',obj.ons);i=i+1;
            str{i}=sprintf('obs infected: %i',obj.oni);i=i+1;
            str{i}=sprintf('obs intensive care: %i',obj.onic);i=i+1;
            str{i}=sprintf('obs vacc: %i',obj.onv);i=i+1;
            str{i}=sprintf('obs dead: %i',obj.ond);i=i+1;
            str{i}=sprintf('obs recovered: %i',obj.onr);i=i+1;
            str{i}=sprintf('obs isolated: %i',obj.oniso);i=i+1;
            
            for k=1:i-1
                disp(str{k})
            end
        end
        %%
        function [ind]=CompileIndizesWithRandomComponent(obj,col,T)
            I=spones(T);
            ri=find(sum(I,2)>1);
            imax=max(max(sum(I,2)));
            
            n=length(ri);
            ind=zeros(imax,n);
            for i=1:n
                res=find(col==ri(i));
                for j=1:length(res)
                    ind(j,i)=res(j);
                end
            end
        end
        
        %%
        function R = CompileRandomMatrix(obj,ti,tj,indr,M,N)
            
            [m,n]=size(indr);
            k=1;
            for i=1:n
                r=rand;
                for j=1:m
                    if(indr(j,i)>0)
                        ri(k)=ti(indr(j,i));
                        rj(k)=tj(indr(j,i));
                        rv(k)=r;
                        k=k+1;
                    end
                end
            end
            
            R=sparse(ri,rj,rv,M,N);
            
        end
        
        %%
        function C=CompileStateChange(obj,T,R,M,N)
            % find the largest negative element row-wise
            
            for i=1:M
                r=max(R(i,:));
                if(~r)
                    l(i)= i;
                    m(i)= find(T(i,:));
                else
                    l(i)= i;
                    ind=find(T(i,:));
                    k=1;
                    r=r-T(i,ind(k));
                    while r>0
                        k=k+1;
                        r=r-T(i,ind(k));
                    end
                    m(i)=ind(k);
                end
            end
            
            C=sparse(l,m,ones(size(l)),M,N);
        end
        
        %%
        function obj=FindAction(obj)
            % greedy algorithm
            Q=zeros(1,obj.na);
            for i=1:obj.na
                Q(i)=obj.R(:,i)'*obj.b;
            end
            ind=find(max(Q)==Q);
            m=length(ind);
            if(m>1)
                r=rand(1,m);
                [~,i]=min(r);
                obj.a=ind(i);
            else
                obj.a=ind;
            end
        end
        %%
        function str=DisplayActions(obj)
            
            str=obj.GetActions;
            
            for k=1:length(str)
                disp(str{k})
            end
        end
        %%
        function str=GetActions(obj)
            i=1;
            str{i}=sprintf('do nothing');i=i+1;
            str{i}=sprintf('detect');i=i+1;
            str{i}=sprintf('vaccinate');i=i+1;
            str{i}=sprintf('isolate infectious');i=i+1;
            str{i}=sprintf('isolate susceptible');i=i+1;
            
            
            if( (i-1) ~= obj.na)
                warning('wrong number of possible action obj.na=%i',obj.na);
            end
        end
        
        %%
        function str=GetParameters(obj)
            
            i=1;
            str{i}=sprintf('disease duration %i days',obj.DiseaseDuration);i=i+1;
            str{i}=sprintf('incubation period %i days',obj.IncubationPeriod);i=i+1;
            str{i}=sprintf('probability of needing intensive care when infected : %.1f%%',100*obj.IntensiveCare);i=i+1;
            str{i}=sprintf('probability of dying in intensive care : %.1f%%',100*(1-obj.IntensiveCareRecovery));i=i+1;
            str{i}=sprintf('days in intensive care : %i ',obj.IntensiveCareTime);i=i+1;
            str{i}=sprintf('probability of not detecting an infectious state wo   test set: %.1f%%',100*obj.notdetect);i=i+1;
            str{i}=sprintf('probability of not detecting an infectious state with test set: %.1f%%',100*obj.notdetectTestSet);i=i+1;
            str{i}=sprintf('probability of sucessful vaccination: %.1f%%',100*obj.vacc);i=i+1;
            str{i}=sprintf('days in isolation : %i ',obj.IsolationDuration);i=i+1;
            str{i}=sprintf('probability of leaving isolation when susceptible: %.1f%%',100*(obj.notisolationSus));i=i+1;
            str{i}=sprintf('probability of leaving isolation when infectious: %.1f%%',100*(obj.notisolationInf));i=i+1;

        end
        %%
        function str=DisplayParameters(obj)
            
            str=obj.GetParameters;
            
            for k=1:length(str)
                disp(str{k})
            end
        end
        
        %%
        function obj=UpdatePerson(obj,p)
            obj.p=p;
            obj.UpdateState;
            obj.UpdateObservation;
            obj.UpdateBelief;
            obj.FindAction;
            obj.UpdateReward;
        end
        
        %%
        function obj=UpdateState(obj)
            v=[];ti=[];tj=[];
            switch(obj.a)
                case 1 % do nothing
                    k=1;
                    % getting infected or not
                    ti(k)=obj.ns;tj(k)=obj.nib;v(k)=obj.p;    k=k+1;
                    ti(k)=obj.ns;tj(k)=obj.ns;v(k)=1-obj.p;  k=k+1;
                    % changing state within the infection
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nib+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of period)
                    ti(k)=obj.nie;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying recovered
                    ti(k)=obj.nr;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying vaccinated
                    ti(k)=obj.nv;tj(k)=obj.nv;v(k)=1;k=k+1;
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nitb+l;
                        v(k)=obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nd;
                        v(k)=1-obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    
                    % recovering (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nr;v(k)=obj.IntensiveCareRecovery;k=k+1;
                    % dying of the infection (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nd;v(k)=1-obj.IntensiveCareRecovery;k=k+1;
                    % isolation infected
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nisob+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of isolation)
                    % or still being infectious
                    d=obj.DiseaseDuration-obj.IsolationDuration;
                    if(d<=0)
                        ti(k)=obj.nisoe;tj(k)=obj.nr;v(k)=1;k=k+1;
                    else
                        ti(k)=obj.nisoe;tj(k)=obj.nie-d;v(k)=1;k=k+1;
                    end
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.nisosusb+l;
                        v(k)=1-obj.notisolationSus;
                        k=k+1;
                    end
                    % leaving isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.ns;
                        v(k)=obj.notisolationSus;
                        k=k+1;
                    end
                    %  back to the susceptible pool
                    ti(k)=obj.nisosuse;tj(k)=obj.ns;v(k)=1;k=k+1;
                    
                    % staying dead
                    ti(k)=obj.nd;tj(k)=obj.nd;v(k)=1;k=k+1;
                    
                    
                    obj.T{obj.a}=sparse(ti,tj,v,obj.n,obj.n);
                    
                    
                    
                case 2 % detect
                    k=1;
                    % getting infected or not
                    ti(k)=obj.ns;tj(k)=obj.nib;v(k)=obj.p;    k=k+1;
                    ti(k)=obj.ns;tj(k)=obj.ns;v(k)=1-obj.p;  k=k+1;
                    % changing state within the infection
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nib+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of period)
                    ti(k)=obj.nie;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying recovered
                    ti(k)=obj.nr;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying vaccinated
                    ti(k)=obj.nv;tj(k)=obj.nv;v(k)=1;k=k+1;
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nitb+l;
                        v(k)=obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nd;
                        v(k)=1-obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    
                    % recovering (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nr;v(k)=obj.IntensiveCareRecovery;k=k+1;
                    % dying of the infection (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nd;v(k)=1-obj.IntensiveCareRecovery;k=k+1;
                    % isolation infected
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nisob+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of isolation)
                    % or still being infectious
                    d=obj.DiseaseDuration-obj.IsolationDuration;
                    if(d<=0)
                        ti(k)=obj.nisoe;tj(k)=obj.nr;v(k)=1;k=k+1;
                    else
                        ti(k)=obj.nisoe;tj(k)=obj.nie-d;v(k)=1;k=k+1;
                    end
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.nisosusb+l;
                        v(k)=1-obj.notisolationSus;
                        k=k+1;
                    end
                    % leaving isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.ns;
                        v(k)=obj.notisolationSus;
                        k=k+1;
                    end
                    %  back to the susceptible pool
                    ti(k)=obj.nisosuse;tj(k)=obj.ns;v(k)=1;k=k+1;
                    % staying dead
                    ti(k)=obj.nd;tj(k)=obj.nd;v(k)=1;k=k+1;
                    
                    
                    obj.T{obj.a}=sparse(ti,tj,v,obj.n,obj.n);
                    
                    
                case 3 % vaccinate
                    k=1;
                    % getting infected or not
                    ti(k)=obj.ns;tj(k)=obj.nib;v(k)=min(obj.p,1-obj.vacc);    k=k+1;
                    ti(k)=obj.ns;tj(k)=obj.ns;v(k)=1-min(obj.p,1-obj.vacc)-min(obj.vacc,1-obj.p);  k=k+1;
                    ti(k)=obj.ns;tj(k)=obj.nv;v(k)=min(obj.vacc,1-obj.p);    k=k+1;
                    % changing state within the infection
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nib+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % becoming recovered (after end of period)
                    ti(k)=obj.nie;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying recovered
                    ti(k)=obj.nr;tj(k)=obj.nr;v(k)=1;k=k+1;
                    
                    
                    
                    % staying vaccinated
                    ti(k)=obj.nv;tj(k)=obj.nv;v(k)=1;k=k+1;
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nitb+l;
                        v(k)=obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nd;
                        v(k)=1-obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    
                    % recovering (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nr;v(k)=obj.IntensiveCareRecovery;k=k+1;
                    % dying of the infection (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nd;v(k)=1-obj.IntensiveCareRecovery;k=k+1;
                    % isolation infected
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nisob+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of isolation)
                    % or still being infectious
                    d=obj.DiseaseDuration-obj.IsolationDuration;
                    if(d<=0)
                        ti(k)=obj.nisoe;tj(k)=obj.nr;v(k)=1;k=k+1;
                    else
                        ti(k)=obj.nisoe;tj(k)=obj.nie-d;v(k)=1;k=k+1;
                    end
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.nisosusb+l;
                        v(k)=1-obj.notisolationSus;
                        k=k+1;
                    end
                    % leaving isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.ns;
                        v(k)=obj.notisolationSus;
                        k=k+1;
                    end
                    %  back to the susceptible pool
                    ti(k)=obj.nisosuse;tj(k)=obj.ns;v(k)=1;k=k+1;
                    % staying dead
                    ti(k)=obj.nd;tj(k)=obj.nd;v(k)=1;k=k+1;
                    
                    
                    obj.T{obj.a}=sparse(ti,tj,v,obj.n,obj.n);
                    
                    % Transition matrix T
                    
                case 4 % isolate infectious
                    k=1;
                    % getting infected or not
                    ti(k)=obj.ns;tj(k)=obj.ns;v(k)=1;    k=k+1;
                    %ti(k)=ns;tj(k)=ns;v(k)=1-p;  k=k+1;
                    % isolate infected
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nisob;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of period)
                    ti(k)=obj.nie;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying recovered
                    ti(k)=obj.nr;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying vaccinated
                    ti(k)=obj.nv;tj(k)=obj.nv;v(k)=1;k=k+1;
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nitb+l;
                        v(k)=obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nd;
                        v(k)=1-obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    
                    % recovering (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nr;v(k)=obj.IntensiveCareRecovery;k=k+1;
                    % dying of the infection (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nd;v(k)=1-obj.IntensiveCareRecovery;k=k+1;
                    
                    % isolation infected
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nisob+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of isolation)
                    % or still being infectious
                    d=obj.DiseaseDuration-obj.IsolationDuration;
                    if(d<=0)
                        ti(k)=obj.nisoe;tj(k)=obj.nr;v(k)=1;k=k+1;
                    else
                        ti(k)=obj.nisoe;tj(k)=obj.nie-d;v(k)=1;k=k+1;
                    end
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.nisosusb+l;
                        v(k)=1-obj.notisolationSus;
                        k=k+1;
                    end
                    % leaving isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.ns;
                        v(k)=obj.notisolationSus;
                        k=k+1;
                    end
                    %  back to the susceptible pool
                    ti(k)=obj.nisosuse;tj(k)=obj.ns;v(k)=1;k=k+1;
                    % staying dead
                    ti(k)=obj.nd;tj(k)=obj.nd;v(k)=1;k=k+1;
                    
                    
                    obj.T{obj.a}=sparse(ti,tj,v,obj.n,obj.n);
                    
                case 5 % isolate susceptible
                    k=1;
                    % getting infected or not
                    ti(k)=obj.ns;tj(k)=obj.nisosusb;v(k)=1;    k=k+1;
                    %ti(k)=ns;tj(k)=ns;v(k)=1-p;  k=k+1;
                    % infected
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nib+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        ti(k)=obj.nib+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of period)
                    ti(k)=obj.nie;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying recovered
                    ti(k)=obj.nr;tj(k)=obj.nr;v(k)=1;k=k+1;
                    % staying vaccinated
                    ti(k)=obj.nv;tj(k)=obj.nv;v(k)=1;k=k+1;
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nitb+l;
                        v(k)=obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    for l=1:obj.IntensiveCareTime
                        ti(k)=obj.nitb+(l-1);
                        tj(k)=obj.nd;
                        v(k)=1-obj.IntensiveCareRecovery;
                        k=k+1;
                    end
                    
                    % recovering (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nr;v(k)=obj.IntensiveCareRecovery;k=k+1;
                    % dying of the infection (after end of period)
                    ti(k)=obj.nite;tj(k)=obj.nd;v(k)=1-obj.IntensiveCareRecovery;k=k+1;
                    % isolation infected
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nisob+l;
                        v(k)=1-obj.IntensiveCare;
                        k=k+1;
                    end
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisob+(l-1);
                        tj(k)=obj.nitb;
                        v(k)=obj.IntensiveCare;
                        k=k+1;
                    end
                    % recovering (after end of isolation)
                    % or still being infectious
                    d=obj.DiseaseDuration-obj.IsolationDuration;
                    if(d<=0)
                        ti(k)=obj.nisoe;tj(k)=obj.nr;v(k)=1;k=k+1;
                    else
                        ti(k)=obj.nisoe;tj(k)=obj.nie-d;v(k)=1;k=k+1;
                    end
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.nisosusb+l;
                        v(k)=1-obj.notisolationSus;
                        k=k+1;
                    end
                    % leaving isolation susceptible
                    for l=1:obj.IsolationDuration
                        ti(k)=obj.nisosusb+(l-1);
                        tj(k)=obj.ns;
                        v(k)=obj.notisolationSus;
                        k=k+1;
                    end
                    %  back to the susceptible pool
                    ti(k)=obj.nisosuse;tj(k)=obj.ns;v(k)=1;k=k+1;
                    % staying dead
                    ti(k)=obj.nd;tj(k)=obj.nd;v(k)=1;k=k+1;
                    
                    
                    obj.T{obj.a}=sparse(ti,tj,v,obj.n,obj.n);
                    
                    
            end
            
            % find transitions with random component
            if(obj.vectorize)
                % compute a realisation for the whole transition matrix,
                % i.e. regardless of the current state only depending on
                % the action
                indr=obj.CompileIndizesWithRandomComponent(ti,obj.T{obj.a});
                R = obj.CompileRandomMatrix(ti,tj,indr,obj.n,obj.n);
                obj.TR = obj.CompileStateChange(obj.T{obj.a},R,obj.n,obj.n);
                
                obj.s=obj.TR'*obj.s;
                
            else
                % compute a realisation only for the row of the transition
                % matrix that corresponds to the actual state
                i=find(obj.s);
                ind=find(obj.T{obj.a}(i,:));
                if(length(ind)>1)
                    j=1;
                    r=rand;
                    r=r-obj.T{obj.a}(i,ind(j));
                    while r>0
                        j=j+1;
                        r=r-obj.T{obj.a}(i,ind(j));
                    end
                    obj.s=zeros(size(obj.s));
                    obj.s(ind(j))=1;
                else
                    obj.s=zeros(size(obj.s));
                    obj.s(ind)=1;
                end
            end
            
            
        end
        %%
        function obj=UpdateObservation(obj)
            v=[];oi=[];oj=[];
            switch(obj.a)
               
                case 2 % detect
                    k=1;
                    %  susceptible
                    oi(k)=obj.ons;oj(k)=obj.ons;v(k)=0.99;    k=k+1;
                    oi(k)=obj.ons;oj(k)=obj.oni;v(k)=0.01;    k=k+1;
                    % detecting an infected
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.oni;
                        v(k)=1-obj.notdetectTestSet;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.ons;
                        v(k)=0.5*obj.notdetectTestSet;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.onr;
                        v(k)=0.5*obj.notdetectTestSet;
                        k=k+1;
                    end
                    oi(k)=obj.nie;oj(k)=obj.oni;v(k)=1-obj.notdetectTestSet;k=k+1;
                    oi(k)=obj.nie;oj(k)=obj.ons;v(k)=0.5*obj.notdetectTestSet;k=k+1;
                    oi(k)=obj.nie;oj(k)=obj.onr;v(k)=0.5*obj.notdetectTestSet;k=k+1;
                    %  recovered
                    oi(k)=obj.nr;oj(k)=obj.onr;v(k)=1.0;k=k+1;
                    
                    % vaccinated
                    oi(k)=obj.nv;oj(k)=obj.onv;v(k)=0.8;k=k+1;
                    oi(k)=obj.nv;oj(k)=obj.onr;v(k)=0.1;k=k+1;
                    oi(k)=obj.nv;oj(k)=obj.ons;v(k)=0.1;k=k+1;
                    
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        oi(k)=obj.nitb+(l-1);
                        oj(k)=obj.onic;
                        v(k)=1;
                        k=k+1;
                    end
                    oi(k)=obj.nite;oj(k)=obj.onic;v(k)=1;k=k+1;
                    % dead
                    oi(k)=obj.nd;oj(k)=obj.ond;v(k)=1;k=k+1;
                    % isolation
                    for l=1:obj.IsolationDuration
                        oi(k)=obj.nisob+(l-1);
                        oj(k)=obj.oniso;
                        v(k)=1;
                        k=k+1;
                    end
                    % last isolation day
                    oi(k)=obj.nisoe;oj(k)=obj.oniso;v(k)=1;k=k+1;
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        oi(k)=obj.nisosusb+(l-1);
                        oj(k)=obj.oniso;
                        v(k)=1;
                        k=k+1;
                    end
                    % last isolation day
                    oi(k)=obj.nisosuse;oj(k)=obj.oniso;v(k)=1;k=k+1;
                    
                    obj.O{obj.a}=sparse(oi,oj,v,obj.n,obj.no);
                    
   
                otherwise 
                    k=1;
                    %  susceptible
                    oi(k)=obj.ons;oj(k)=obj.ons;v(k)=0.98;    k=k+1;
                    oi(k)=obj.ons;oj(k)=obj.oni;v(k)=0.01;    k=k+1;
                    oi(k)=obj.ons;oj(k)=obj.oniso;v(k)=0.01;    k=k+1;
                    % detecting an infected
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.oni;
                        v(k)=1-obj.notdetect;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.ons;
                        v(k)=0.5*obj.notdetect;
                        k=k+1;
                    end
                    for l=1:obj.DiseaseDuration
                        oi(k)=obj.nib+(l-1);
                        oj(k)=obj.onr;
                        v(k)=0.5*obj.notdetect;
                        k=k+1;
                    end
                    oi(k)=obj.nie;oj(k)=obj.oni;v(k)=1-obj.notdetect;k=k+1;
                    oi(k)=obj.nie;oj(k)=obj.ons;v(k)=0.5*obj.notdetect;k=k+1;
                    oi(k)=obj.nie;oj(k)=obj.onr;v(k)=0.5*obj.notdetect;k=k+1;
                    %  recovered
                    oi(k)=obj.nr;oj(k)=obj.onr;v(k)=1.0;k=k+1;
                    % vaccinated
                    oi(k)=obj.nv;oj(k)=obj.onv;v(k)=0.8;k=k+1;
                    oi(k)=obj.nv;oj(k)=obj.onr;v(k)=0.1;k=k+1;
                    oi(k)=obj.nv;oj(k)=obj.ons;v(k)=0.1;k=k+1;
                    
                    % intensive care
                    for l=1:obj.IntensiveCareTime
                        oi(k)=obj.nitb+(l-1);
                        oj(k)=obj.onic;
                        v(k)=1.0;
                        k=k+1;
                    end
                    oi(k)=obj.nite;oj(k)=obj.onic;v(k)=1.0;k=k+1;
                    
                    %
                    oi(k)=obj.nd;oj(k)=obj.ond;v(k)=1;k=k+1;
                    % isolation
                    for l=1:obj.IsolationDuration
                        oi(k)=obj.nisob+(l-1);
                        oj(k)=obj.oniso;
                        v(k)=1;
                        k=k+1;
                    end
                    % last isolation day
                    oi(k)=obj.nisoe;oj(k)=obj.oniso;v(k)=1;k=k+1;
                    % isolation susceptible
                    for l=1:obj.IsolationDuration
                        oi(k)=obj.nisosusb+(l-1);
                        oj(k)=obj.oniso;
                        v(k)=1;
                        k=k+1;
                    end
                    % last isolation day
                    oi(k)=obj.nisosuse;oj(k)=obj.oniso;v(k)=1;k=k+1;
                    
                    obj.O{obj.a}=sparse(oi,oj,v,obj.n,obj.no);
                    
            end
            
            if(obj.vectorize)
                % same as for transition matrix
                oindr=obj.CompileIndizesWithRandomComponent(oi,obj.O{obj.a});
                oR = obj.CompileRandomMatrix(oi,oj,oindr,obj.n,obj.no);
                obj.OR = obj.CompileStateChange(obj.O{obj.a},oR,obj.n,obj.no);
                obj.z=obj.OR'*obj.s;
            else
                % compute a realisation only for the row of the transition
                % matrix that corresponds to the actual state
                i=find(obj.s);
                ind=find(obj.O{obj.a}(i,:));
                if(length(ind)>1)
                    j=1;
                    r=rand;
                    r=r-obj.O{obj.a}(i,ind(j));
                    while r>0
                        j=j+1;
                        r=r-obj.O{obj.a}(i,ind(j));
                    end
                    obj.z=zeros(size(obj.z));
                    obj.z(ind(j))=1;
                else
                    obj.z=zeros(size(obj.z));
                    obj.z(ind)=1;
                end
                
            end
            
        end
        %%
        function obj=UpdateBelief(obj)
            W=obj.O{obj.a}(:,find(obj.z)).*obj.T{obj.a}';
            bs=W*obj.b;
            
            sbs=sum(bs);
            if(sbs>0)
                bs=bs/sbs;
            else
                warning('Something went HORROBLY wrong s=%i,o=%i,a=%i. Please check transision and observation matrizes for consistency. Continuing with flat belief.',find(obj.s),find(obj.z),obj.a)
                bs=ones(size(bs));
                bs=bs/sum(bs);
            end
            obj.b=bs;
        end
        
        
        function n=GetNumberOfStates(obj)
            n=obj.n;
        end
        
        function n=GetNumberOfObservations(obj)
            n=obj.no;
        end
        
        function n=GetNumberOfActions(obj)
            n=obj.na;
        end
        
        
        function status=IsSusceptible(obj)
            if(obj.s(obj.ns)==1)
                status=true;
            else
                status=false;
            end
        end
        
        %%
        function status=SpreadsInfection(obj)
            i=find(obj.s);
            if( (i>=obj.nib && i<=obj.nie) )
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsInfectious(obj)
            i=find(obj.s);
            if( (i>=obj.nisob && i<=obj.nisoe) || (i>=obj.nib && i<=obj.nie) )
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsRecovered(obj)
            if(obj.s(obj.nr)==1)
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsVaccinated(obj)
            if(obj.s(obj.nv)==1)
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsDead(obj)
            if(obj.s(obj.nd)==1)
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsInIntensiveCare(obj)
            i=find(obj.s);
            if(i>=obj.nitb && i<=obj.nite)
                status=true;
            else
                status=false;
            end
        end
        %%
        function status=IsInIsolation(obj)
            i=find(obj.s);
            if( (i>=obj.nisob && i<=obj.nisoe) ||  (i>=obj.nisosusb && i<=obj.nisosuse))
                status=true;
            else
                status=false;
            end
        end
        %%
        function obj=InitRewardMatrix(obj)
            % R(state,action)
            obj.R=zeros(obj.n,obj.na);
            if(~isempty(obj.DNA))
                % load a DNA sequence into the reward matrix
                if(length(obj.DNA)==obj.na*obj.n)
                    for i=1:obj.n
                        for j=1:obj.na
                            obj.R(i,j)=obj.DNA((i-1)*obj.na+j);
                        end
                    end
                else
                    error('wrong genom size in reward matrix initialization');
                end
            else
                if(true)
                    %1    'do nothing'
                    %2    'detect'
                    %3    'vaccinate'
                    %4    'isolate infectious'
                    %5    'isolate susceptible'
                    
                    % reward should be b(p), i.e. belief of probability
                    % of getting infected dependent.
                    % vaccination should be obj.vacc dependent
                    
                    % susceptible
                    obj.R(obj.ns,1)=0;
                    obj.R(obj.ns,2)=0;
                    obj.R(obj.ns,3)=obj.vacc;
                    obj.R(obj.ns,4)=0;
                    obj.R(obj.ns,5)=0;
                    % infectious
                    obj.R(obj.nib:obj.nie,1)=0;
                    obj.R(obj.nib:obj.nie,2)=2;
                    obj.R(obj.nib:obj.nie,3)=0;
                    obj.R(obj.nib:obj.nie,4)=5;
                    obj.R(obj.nib:obj.nie,5)=0;
                    % recovered
                    obj.R(obj.nr,1)=1;
                    obj.R(obj.nr,2)=0;
                    obj.R(obj.nr,3)=0;
                    obj.R(obj.nr,4)=0;
                    obj.R(obj.nr,5)=0;
                    % vaccinated
                    obj.R(obj.nv,1)=1;
                    obj.R(obj.nv,2)=0;
                    obj.R(obj.nv,3)=0;
                    obj.R(obj.nv,4)=0;
                    obj.R(obj.nv,5)=0;
                    % intensive care
                    obj.R(obj.nitb:obj.nite,1)=1;
                    obj.R(obj.nitb:obj.nite,2)=0;
                    obj.R(obj.nitb:obj.nite,3)=0;
                    obj.R(obj.nitb:obj.nite,4)=0;
                    obj.R(obj.nitb:obj.nite,5)=0;
                    % isolated infectious
                    obj.R(obj.nisob:obj.nisoe,1)=1;
                    obj.R(obj.nisob:obj.nisoe,2)=0;
                    obj.R(obj.nisob:obj.nisoe,3)=0;
                    obj.R(obj.nisob:obj.nisoe,4)=0;
                    obj.R(obj.nisob:obj.nisoe,5)=0;
                    % isolated susceptible
                    obj.R(obj.nisosusb:obj.nisosuse,1)=1-obj.vacc;
                    obj.R(obj.nisosusb:obj.nisosuse,2)=0;
                    obj.R(obj.nisosusb:obj.nisosuse,3)=obj.vacc;
                    obj.R(obj.nisosusb:obj.nisosuse,4)=0;
                    obj.R(obj.nisosusb:obj.nisosuse,5)=0;
                    % dead is very bad
                    obj.R(obj.nd,1)=1;
                    obj.R(obj.nd,2)=0;
                    obj.R(obj.nd,3)=0;
                    obj.R(obj.nd,4)=0;
                    obj.R(obj.nd,5)=0;
                    
                    %obj.R(:,:)=0;
                end
            end
        end
        
        %%
        function obj=DisplayRewardMatrix(obj)
            figure
            clf
            str=obj.DisplayActions;
            for i=1:obj.na
                subplot(obj.na,1,i)
                plot(obj.R(:,i))
                title(str{i})
                xlabel('state')
                ylabel('reward')
            end
        end
        %%
        function obj=UpdateReward(obj)
            obj.reward=obj.reward+obj.R(find(obj.s),obj.a);
        end
        %%
        function obj=ConsistencyCheck(obj)
            a_save=obj.a;
            for i=1:obj.na
                obj.a=i;
                obj.UpdateState;
                tp=sum(obj.T{obj.a},2);
                if(find(tp~=1))
                    figure
                    clf
                    subplot(2,1,1)
                    plot(tp)
                    subplot(2,1,2)
                    spy(obj.T{obj.a})
                    error('transition probabilities do not add up to 1 (a=%i)',obj.a)
                end
                obj.UpdateObservation;
                op=sum(obj.O{obj.a},2);
                if(find(op~=1))
                    figure
                    clf
                    subplot(2,1,1)
                    plot(op)
                    subplot(2,1,2)
                    spy(obj.O{obj.a})
                    error('observation probabilities do not add up to 1 (a=%i)',obj.a)
                end
            end
            obj.a=a_save;
        end
        
    end
end

