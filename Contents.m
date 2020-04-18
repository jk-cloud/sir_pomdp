% SIR POMDP toolbox  3/2020
%
% toolbox for simulation the spread of a disease and choosing appropriate
% actions to hinder the spread. 
%
% The disease is modeled in compartments
%   Susceptible
%   Infectioussss
%   Recovered
%
% In order to model actions as well compartments like
%   vaccinated
%   isolated being infectious
%   isolated being susceptible
%   intensive care
%   dead 
%   ...
%  are added. The transitions of changes between these
%  compartments, i.e. states, are modeled using transition matrices.
%
%  Possible observations are also modeled.
%  
%  The model uses a POMDP approach to model the uncertainty in
%  observations (Partially Observable Markov Decision Process).
%  I recommend to read the article "planning and acting in partially
%  observable stochastic domains" by L.P. Kaebling et al. 1997.
%
%  The possible actions are 
%   do nothing
%   try to find out if someone is infectious
%   isolate an infectious
%   isolate a susceptible
%   vaccinate
%
%   Given a set of "rewards" for each state and action combination
%   the model tries to find the optimal "greedy" action.
%

%
%  (c) 2020 Jens Kappey and the sir_pomdp contributors.
%