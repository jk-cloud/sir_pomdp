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

status = false;


P=Person(0.1,DNA);
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
    plot(a)
    xlabel('time')
    ylabel('action')
    
    
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
    
    
    P.DisplayRewardMatrix
end
status = true;
end

