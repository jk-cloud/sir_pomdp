function  [S,I,R,V,ISO,D,ICA,REWARD] = test_stoc_age_community(DNA)
%TEST_STOC_AGE_COMMUNITY test routine for stochastic simulation of a community
%   Driver routine for simulating the spread of a disease within a 
%   community of people using many realistations of the function
%
%    stoc_community
% 
%   to build averages.
%
%  usage:  [S,I,R,V,ISO,D,ICA,REWARD] = test_stoc_age_community([DNA])
%
%      where DNA sets the reward matrix in the class Person
%      DNA con be used to interface the code using a genetic algorithm.
%      in order to optimze the reward matrix.s

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%
switch(nargin)
    case 1
      
      PLOT=true;
    otherwise
      DNA=[];
      PLOT=true;
end

M=6;      % realisations  
N=10000;   % PopulationSize
steps=30; % time steps

NeedForIntensiveCareAge=[1,20,50,100];
NeedForIntensiveCare   =[0.001,0.001,0.01,0.1];
IntensiveCareRecoveryAge=[1,20,50,100];
IntensiveCareRecoveryProb=[0.99,0.95,0.9,0.6];
IntensiveCareTimeAge=[1,20,35,50,100];
IntensiveCareTimeDays=[10,11,12,15,25];

p=Person(0);
p.ConsistencyCheck;

%DNA=3*rand(1,p.GetNumberOfActions*p.GetNumberOfStates);
%DNA=round(DNA);

C = community(DNA);
C.SetPopulationSize(N);
C.SetSimulationSteps(steps);
C.SetNumberOfRealizations(M);
C.SetNeedForIntensiveCareBasedOnAge(NeedForIntensiveCareAge,NeedForIntensiveCare);
C.SetIntensiveCareRecoveryBasedOnAge(IntensiveCareRecoveryAge,IntensiveCareRecoveryProb);
C.SetIntensiveCareTimeBasedOnAge(IntensiveCareTimeAge,IntensiveCareTimeDays);
C.Initialize;
C.Evolve;

[S,I,R,V,ISO,D,ICA,REWARD,A] = C.ReturnResults;



if(PLOT)
    phi=-20;
    psi=30;
    figure(4)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','evolution of a community (averaged)')
    clf
    subplot(9,1,1)
     set(gca,'Visible','off')
    txt=p.GetParameters;
    l=length(txt);
    n=ceil(l/2);
    txt1=txt(1:n);
    txt2=txt(n+1:l);
    text(0.1,0.5,txt1,'FontSize',14)
    text(0.6,0.5,txt2,'FontSize',14)
    
    subplot(9,1,2)
    hold on
    plot(S','.')
    pl=plot(mean(S),'r-');
    set(pl,'LineWidth',3);
    ylabel('susceptible')
    
    
    subplot(9,1,3)
    hold on
    plot(I','.')
    pl=plot(mean(I),'r-');
    set(pl,'LineWidth',3);
    ylabel('infectious')
    
    subplot(9,1,4)
    hold on
    plot(R','.')
    pl=plot(mean(R),'r-');
    set(pl,'LineWidth',3);
    ylabel('recovered')
    
    subplot(9,1,5)
    hold on
    plot(V','.')
    pl=plot(mean(V),'r-');
    set(pl,'LineWidth',3);
    ylabel('vaccinated')
    
    subplot(9,1,6)
    hold on
    plot(ISO','.')
    pl=plot(mean(ISO),'r-');
    set(pl,'LineWidth',3);
    
    ylabel('isolated')
    
    subplot(9,1,7)
    hold on
    plot(D','.')
    pl=plot(mean(D),'r-');
    set(pl,'LineWidth',3);
    
    ylabel('dead')
    
    subplot(9,1,8)
    hold on
    plot(ICA','.')
    pl=plot(mean(ICA),'r-');
    set(pl,'LineWidth',3);
    ltxt=sprintf('∫ = %.1f ',trapz(mean(ICA)));
    ylabel('intensive care')
    legend(ltxt)
    
    subplot(9,1,9)
    hold on
    plot(REWARD','.')
    pl=plot(mean(REWARD),'r-');
    set(pl,'LineWidth',3);
    xlabel('time [day]')
    ylabel('total reward')
    
end
end




