# TraitDistributionEvolution
Code relating to paper titled "Distinguishing Passive and Driven Trait Evolution in the presence of Boundaries".

MOMv3.3.txt contains all of the mammal mass data analysed in this paper.
The analysis of these data is performed by the code in ProcessingMammalData.m.

BivalveSummary.csv and BrachiopodSummary.csv contain all of the bivalve and brachiopod linear size data analysed in this paper.
The analysis of these data is performed by the code in BrachiopodsBivalvesAnalysisMatched.m.

The remaining functions are required to run all of the analysis in these two data analysis codes.

SimulationDis.m is the exception - this is the discrete particle based simulation which is based on the following dynamics: particles (representing species) can undergo speciation (particle birth), extinction (particle death), and evolution (particle trait change).  The evolution can occur via either a symmetric mechanism (the probability is equal to increase or decrease the trait value of the particle) or a driven mechanism (the particle's trait changes with 100% certainty in a certain direction - dictated by the sign of the parameter a).  The simulation is done in discrete trait space and so particles are confined to a discrete number of bins along the trait axis.  The simulation is performed via a Gillespie style algorithm.  There is an option to run the simulation multiple times and average the results.  There is also the option to change the rules at the boundaries on the trait space.  The possible boundary types are a Dirichlet boundary, a Neumann boundary and a zero-flux (Robin-type) boundary.  
