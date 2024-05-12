# Experience dependent behavior in C. elegans

The simulations of Figure 5, 6 and 7 of the paper are generated from figures_paper_step.ipynb.
The simulations of Figure 9, 10 and 11 of the paper are generated from graphs_paper.ipynb.

To run the simulation in C++
1. run 
mingw32-make
2. run
./main
3. for cleaning 
mingw32-make clean

In the main.cpp file you can choose running example trajectories (#define trajectories). In this case, you will save a txt file with the column values detailed in the function PrintDetail in WormAgent.cpp for every time step of the simulation. The second option is to run several simulations (#define chem_index) and save only the last time point of the simulation to compute the chemotaxis index afterwards (with the chemotaxis_index.ipynb notebook). In this case, also a txt file is generated from the PrintDetail function.
 
The NaCl concentration during cultivation can be changed with the parameter "phenotype(25)" (C_breed). We explored the values C_breed=25,50,100.

WormAgent.h
You can choose between a step of NaCl (#define GRAD_STEP), or the gradient use in the paper (#define GRAD_GAUS). Beware that if you use the step concentration, you need to change the constant Preexposure to a value lower than 5000 for memory reasons.



