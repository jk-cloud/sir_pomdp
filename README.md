# sir_pomdp
Matlab Toolbox for disease control using partially observable decision 
process framework for action selection

## Why?

During the Corona pandemic in 2020 with no vaccine available several 
measures where proposed to deal with the epidemic.
 
Which of these measures/actions or combination thereof is most efficient?

Having worked in epidemilogy decades ago I knew how to model a disease.
Having worked with partially obversavble Markov chains also decades ago
I also knew how to find a best set of actions, given that the state of the
world is not well known. So I decided to combine both to find an answer.s

## How?
 The disease is modeled in compartments
  - Susceptible
  - Infectious
  - Recovered

 In order to model actions as well compartments like
  - vaccinated
  - isolated being infectious
  - isolated beiing susceptible
  - intensive care
  - dead 
  - ...
  
  are added. The transistions of changes between these
  compartments, i.e. states, are modelled using transition matrizes.

  Possible observations are also modelled.
  
  The model uses a POMDP approach to model the uncertainty in
  observations (Partially Observable Markov Decision Process).
  I recommend to read the article by "planning and acting in partially
  observabl stochastic domains" by L.P. Kaebling et al. 1997.

  The possible actions are 
  - do nothing
  - try to find out if someone is infectious
  - isolate an infectious
  - isolate a susceptible
  - vaccinate

   Given a set of "rewards" for each state and action combination
   the model tries to find the optimal "greedy" action.
