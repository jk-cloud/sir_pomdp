function  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community(DNA)
%TEST_COMMUNITY test routine for  community
%
%   usage:  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community([DNA])

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%

switch(nargin)
    case 1
      N=1000;    % PopulationSize
      steps=60;
      PLOT=true;

    otherwise
      N=1000;    % PopulationSize
      steps=30;     % time steps  
      DNA=[];
      PLOT=true;
end

NeedForIntensiveCareAge=[1,20,50,100];
NeedForIntensiveCare   =[0.001,0.001,0.01,0.1];
IntensiveCareRecoveryAge=[1,20,50,100];
IntensiveCareRecoveryProb=[0.99,0.95,0.9,0.6];
IntensiveCareTimeAge=[1,20,35,50,100];
IntensiveCareTimeDays=[10,11,12,15,25];



P0=Person(0);
P0.ConsistencyCheck;

C = community(DNA);
C.SetPopulationSize(N);
C.SetSimulationSteps(steps);
C.SetNeedForIntensiveCareBasedOnAge(NeedForIntensiveCareAge,NeedForIntensiveCare);
C.SetIntensiveCareRecoveryBasedOnAge(IntensiveCareRecoveryAge,IntensiveCareRecoveryProb);
C.SetIntensiveCareTimeBasedOnAge(IntensiveCareTimeAge,IntensiveCareTimeDays);
C.Initialize;
C.Evolve;

[S,I,R,V,ISO,D,ICA,REWARD,A] = C.ReturnResults;


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
    set(gca,'Visible','off')
    
    txt=P0.GetParameters;
    l=length(txt);
    n=ceil(l/2);
    txt1=txt(1:n);
    txt2=txt(n+1:l);
    text(0.1,0.5,txt1,'FontSize',14)
    text(0.6,0.5,txt2,'FontSize',14)
    
    
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
    ltxt=sprintf('∫ = %.1f ',trapz(ICA));
    ylabel('intensive care')
    legend(ltxt)
    
    subplot(9,1,9)
    plot(REWARD)
    xlabel('time [day]')
    ylabel('total reward')
    
end
end




