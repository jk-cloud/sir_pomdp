function status = test_person(DNA)
%TEST_PERSON testing the person class
%  
%   usage:  test_person([DNA])

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%


p=Person(0.0);
p.ConsistencyCheck;

switch(nargin)
    case 1
        
    otherwise    
     DNA=3*rand(1,p.GetNumberOfActions*p.GetNumberOfStates);
     DNA=round(DNA);           
end
DNA=[];
status = false;

age= round(100*rand(1));
P=Person(0.1,DNA,age);
% set age variable parameters
NeedForIntensiveCareAge=[1,20,50,100];
NeedForIntensiveCare   =[0.001,0.001,0.01,0.1];
P.SetNeedForIntensiveCareBasedOnAge(NeedForIntensiveCareAge,NeedForIntensiveCare);
IntensiveCareRecoveryAge=[1,20,50,100];
IntensiveCareRecoveryProb=[0.99,0.95,0.9,0.6];
P.SetIntensiveCareRecoveryBasedOnAge(IntensiveCareRecoveryAge,IntensiveCareRecoveryProb);
IntensiveCareTimeAge=[1,20,35,50,100];
IntensiveCareTimeDays=[10,11,12,15,25];
P.SetIntensiveCareTimeBasedOnAge(IntensiveCareTimeAge,IntensiveCareTimeDays);
P.SetParametersBasedOnAge;

P.DisplayTransitionIndizes;
P.DisplayObservationIndizes;

N=30;
S=zeros(P.GetNumberOfStates,N);
B=S;
Z=zeros(P.GetNumberOfObservations,N);
a=zeros(1,N);

for k=1:N
    P.UpdatePerson(0.2);
    Z(:,k)=P.z;
    S(:,k)=P.s;    
    B(:,k)=P.b;
    a(k)  =P.a;
end

PLOT=true;
if(PLOT)
    phi=-20;
    psi=30;
    figure(3)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','evolution of a single person')
    clf
    subplot(1,4,1)
    plot(a,'o')
    xlabel('time')
    ylabel('action')
    lbs=P.GetActions;
    yt=linspace(1,length(lbs),length(lbs));
    yticks(yt);
    yticklabels(lbs);
    
    subplot(1,4,2)
    surf(S)
    xlabel('time')
    ylabel('state')
    view(phi,psi);
    title('evolution')
    
    subplot(1,4,3)
    surf(B-S)
    xlabel('time')
    ylabel('belief-state')
    view(phi,psi);
    title('evolution')
    
    subplot(1,4,4)
    
    surf(B)
    xlabel('time')
    ylabel('belief')
    view(phi,psi);
    title('evolution')
    
    figure(4)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','age dependend properpties of a single person')
    clf
    subplot(3,1,1)
    P.PlotNeedForIntensiveCareBasedOnAge;
    subplot(3,1,2)
    P.PlotIntensiveCareRecoveryBasedOnAge;
    subplot(3,1,3)
    P.PlotIntensiveCareTimeBasedOnAge;
    %% will be displayed in a "next" window
    P.DisplayRewardMatrix
end
status = true;
end

