# MC_LatticeGas
This project is done in 2014, from the Undergraduate Research Opportunities Program (UROP) that I participated during my undergraduate study in HKUST. It is a simplified Monte Carlo (MC) simulation of Br atoms on Cu(111) surface based on lattice gas model (2-D), with 2 competing forces: bonding force between the nearest neighbor and repulsive force (∝1/r^3) between each of atom, without time-related factors like transition rates, thus it only searches the most energetically favorable state for the configuration of atoms under two competing forces.  
The results are stunning. 3 different coverages generate different patterns. The most interesting case is when 33.65% of the surface is covered, the deposition of Br on Cu surface shows [turing patterns](https://en.wikipedia.org/wiki/Turing_pattern) with characteristic stripes in experiment, while my simulation shows the same turing patterns explained by [reaction-diffusion theory](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system). The local chemical reactions are compete with the global diffusion just like the local bonding is compete with global repulsive force though there is no chemical reaction in my system. 
