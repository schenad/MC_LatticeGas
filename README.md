# MC_LatticeGas
This project is from the Undergraduate Research Opportunities Program (UROP) that I participated during my undergraduate study in HKUST (2013-2014). It is a simplified Monte Carlo (MC) simulation of Br atoms on Cu(111) surface based on [lattice gas model](https://en.wikipedia.org/wiki/Ising_model#Lattice_gas) (2-D), with 2 competing forces: bonding force between the nearest neighbor and repulsive force between each of atom (any distance). [Periodic boundary condition](https://en.wikipedia.org/wiki/Periodic_boundary_conditions) is used to eliminate edge effect. This toy model only searches the most energetically favorable state for the configuration of atoms under two competing forces. Searching for the best parameters for the model was done by calculating the size distribution, shape distribution and pair correlation function of the islands in simulated model and comparing these characteristics with the ones in actual experiments. 

## Simulation vs. Experiment

Below shows selected simulation results in compare to the experiment:

- 10.92% Case in Experiment (STM): large clusters are molecules that are irrelevant to our Br/Cu system.

![Experiment of Low Coverage Case](/Results/Example_STM_Scans/10%25_Z%20TraceUp%20Tue%20May%2028%2015_48_39%202013%20%5B20-1%5D%20%20STM_AtomManipulation%20STM.bmp)

- 10.92% Case in Simulation with the Optimized Parameters: 

![Simulation of Low Coverage Case](/Results/GIF/10.92%25%2C%20E%3D2e-3%2CA%3D4e-4.gif)

- 33.65% Case in Experiment (STM):

![Experiment of Mediun Coverage Case](/Results/Example_STM_Scans/33%25_Z%20TraceUp%20Tue%20May%2028%2019_23_23%202013%20%5B50-1%5D%20%20STM_AtomManipulation%20STM.bmp)

- 33.65% Case in Simulation with the Optimized Parameters:

![Simulation of Medium Coverage Case](/Results/GIF/33.65%25%2C%20E%3D2e-3%2CA%3D4e-4.gif)

The results are stunning. 3 different coverages generate 2 different patterns (the 10.92% is very similar to 71.95%). The most interesting case is when 33.65% of the surface is covered, the deposition of Br on Cu surface shows [turing patterns](https://en.wikipedia.org/wiki/Turing_pattern) with characteristic stripes in experiment, while my simulation shows very similar patterns. The turing pattern can be explained by [reaction-diffusion theory](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system). The local chemical reactions are compete with the global diffusion just like the local bonding is compete with global repulsive force. Although there is no chemical reaction in my system, similar pattern is formed and serves as a great example of [pattern formation](https://en.wikipedia.org/wiki/Pattern_formation) with simple rules. In fact, in [classical nucleation theory](https://en.wikipedia.org/wiki/Classical_nucleation_theory), the density of nucleated island mainly depends on the competition between deposition flux (F) and diffusion energy (Ed), which may also explain the formation of turing pattern in my simulation.  

## Phase Transition

With the model in hand, I calculate the critical temperature (Tc) of phase transition for our system. I calculate for the cases with and without the repulsion for different coverages. Without repulsion, the phase transition is significant at Tc as shown below, while the repulsive force gives a more gradual transition. The specific heat capacity (Cv) as a funtion of temperature was calculated to characterize the phase transistion. 

- Coverage = 1%:

![1%](/Results/Cv%20Calculation/No%20Repulsion/1%25.png)

- Coverage = 5%:

![5%](/Results/Cv%20Calculation/No%20Repulsion/5%25.png)

- Coverage = 10.92%:

![10.92%](/Results/Cv%20Calculation/No%20Repulsion/10.92%25.png)

- Coverage = 33.65%:

![33.65%](/Results/Cv%20Calculation/No%20Repulsion/33.65%25.png)

- Coverage = 50%:

![50%](/Results/Cv%20Calculation/No%20Repulsion/50%25.png)

- Coverage = 71.95%:

![71.95%](/Results/Cv%20Calculation/No%20Repulsion/71.95%25.png)

- Coverage = 90%:

![90%](/Results/Cv%20Calculation/No%20Repulsion/90%25.png)

- Coverage = 99%:

![99%](/Results/Cv%20Calculation/No%20Repulsion/99%25.png)

## Publication

The whole work has been published as a poster presentation (**[click to view](/HKPS%202014.pdf)**) in 17th Annual Conference of the Physical Society of Hong Kong, 06/07/2014 and received the Best Student Poster Award on that year.
