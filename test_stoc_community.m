function  [S,I,R,V,ISO,D,ICA,REWARD] = test_stoc_community(DNA)
%TEST_STOC_COMMUNITY test routine for stochastic simulation of a community
%   Driver routine for simulating the spread of a disease within a 
%   community of people using many realistations of the function
%
%    stoc_community
% 
%   to build averages.
%
%  usage:  [S,I,R,V,ISO,D,ICA,REWARD] = test_stoc_community([DNA])
%
%      where DNA sets the reward matrix in the class Person
%      DNA con be used to interface the code using a genetic algorithm.
%      in order to optimze the reward matrix.s

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%
switch(nargin)
    case 1
      
      PLOT=false;
    otherwise
      DNA=[];
      PLOT=true;
end


p=Person(0);
p.ConsistencyCheck;

%DNA=3*rand(1,p.GetNumberOfActions*p.GetNumberOfStates);
%DNA=round(DNA);

[S,I,R,V,ISO,D,ICA,REWARD] = stoc_community(DNA);

if(PLOT)
    phi=-20;
    psi=30;
    figure(4)
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','evolution of a community (averaged)')
    clf
    subplot(9,1,1)
    
    
    
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





