# SIR_Stochastic

A simple [Susceptable - Infected - Recovered](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology) (SIR) [stochatsic model](https://en.wikipedia.org/wiki/Gillespie_algorithm#Another_example:_The_SIR_epidemic_without_vital_dynamics) for epidemic modelling.

This method has sometimes been called [Dynmamic Monte Carlo](https://en.wikipedia.org/wiki/Dynamic_Monte_Carlo_method), and is 
used in reaction chemistry for predicting the population of some chemical
compounds at a future time, for some copuled reaction pathways with some rate
at which the reactions occur; and that the chemicals are well mixed.

The problem may also be solve with coupled ordinary differential equations, 
but becomes difficult if the ODE's are too stiff. Dynamic MC does not have
this problem.

The problem is almost idential for an SIR model, and works in the following way.
The initial populations of N_I, N_S are given, with N_R = 0; and also
the rates at which people are infected and recovered (InfRate, CureRate).
Starting from t=0, either an infection or a recovery will occor. The time to 
one of these reactions is sampled from an exponential distribution, whose parameter
is calculated from N_I, N_S, InfRate and CureRate. The time is stepped by the sample
and then one of the reactions is sampled. The populations is then ajusted according
to which reaction has occured. This is repeated until t = Tfinal or until N_I == 0.

Simple!           

Because this is a stochastic simulation, the entire simulation may be re-run (batched)
a number of times calculating the variation due to the stochasticity.
In this model this is the Nbatches parameter.

We can also model social distancing. This is parameter V, and it reduces
the rate of infection.

If you find that this model is too slow, we can run it in parallel, 
implement importance sampling, or both.

To try out the stochastic model: runBatchesSIR.m

Uncertainty Propagation
---

Uncertainty (probability distribution) in the infection rate, recovery rate and the spacial parameter may be propagated with Monte Carlo. The Input distributions are Gaussian, but may be anything. Dependency is modelled using a gaussian [copula](https://en.wikipedia.org/wiki/Copula_(probability_theory)), which is sampled with Cholesky decomposition. You may define the correlations in terms of the partial correlations, which
may take any value in [-1, 1], independently.

Because we have both aleatory (stochastic model) and epistemic uncertainty (uknown rates),
the output of this process will be a 2nd order distribution, a distribution of
distributions. We will have one for every point in time. You may take the bounds of the 2nd
order distribution and create a [p-box](https://en.wikipedia.org/wiki/Probability_box). 

Running "runSIRUQ.m" will:

![alt text](https://github.com/AnderGray/SIR_Stochastic/blob/master/Plots/Process.png "All samples of the 2nd order distribution of the process")

Where red, blue and green are the infected, recovered and susceptable numbers respectively.

All of the samples of the 2nd order distribution have been projected onto the same axis. Since a 2nd order distribution is produced, the mean will also have a distribution. The black lines is the "mean of the mean".

A slice of the 2nd order distribution at a specific time may be plotted using "sliceTime.m"

![alt text](https://github.com/AnderGray/SIR_Stochastic/blob/master/Plots/Infected0.7.png "2nd order distribution of Infected numbers at T= 0.7")

![alt text](https://github.com/AnderGray/SIR_Stochastic/blob/master/Plots/Suceptable0.7.png "2nd order distribution of Suceptable numbers at T= 0.7")

![alt text](https://github.com/AnderGray/SIR_Stochastic/blob/master/Plots/Recovered0.7.png "2nd order distribution of Recovered numbers at T= 0.7")


Where the black lines are now the "mean distribution".

You may also produce a time-lapse using "rollingPbox.m"

![alt text](https://github.com/AnderGray/SIR_Stochastic/blob/master/Plots/RollingPboxQuick.gif "2nd order distribution time-lapse of Infected numbers")


Stay home and keep coding!
---

