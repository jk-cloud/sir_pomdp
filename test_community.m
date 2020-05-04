function  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community(N,steps,DNA)
%TEST_COMMUNITY test routine for  community
%
%   usage:  [S,I,R,V,ISO,D,ICA,REWARD,A] = test_community(N,steps,[DNA])

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%

switch(nargin)
    case 1
      steps=30;
      DNA=[];
      PLOT=true;
    case 2
      DNA=[];
      PLOT=true;
    case 3
      PLOT=true;

    otherwise
      N=1000;    % PopulationSize
      steps=30;     % time steps  
      DNA=[];
      PLOT=true;
end

P0=Person(0);
P0.ConsistencyCheck;

[S,I,R,V,ISO,D,ICA,REWARD,A] = community(N,steps,DNA);

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
    ltxt=sprintf('∫ = %.1f ',trapz(ICA));
    ylabel('intensive care')
    legend(ltxt)
    
    subplot(9,1,9)
    plot(REWARD)
    xlabel('time [day]')
    ylabel('total reward')
    
end
end




