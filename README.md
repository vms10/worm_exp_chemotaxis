# Experience dependent behavior in C. elegans
 Code to run the simulations proposed in: María Sol Vidal-Saez, Oscar Vilarroya, and Jordi Garcia-Ojalvo. A Multiscale Sensorimotor Model of Experience-Dependent Behavior in a Minimal Organism.

## Organization

 - The simulations of Figure 5, 6 and 7 of the paper are generated from figures_paper_step.ipynb.
 - The simulations of Figure 9, 10 and 11 of the paper are generated from graphs_paper.ipynb.

## Simulations in C++
The chemotaxis assay simulations run in C++. To run them:

1. Compile the files (with make)

```
mingw32-make
```

2. Run the executable that was generated in the previous step

```
./main
```

3. If you want to delete the .o and .exe files and re-compile, run:

```
mingw32-make clean
```

## Other considerations

### main.cpp
In the main.cpp file you can choose running example trajectories (#define trajectories). In this case, you will save a txt file with the column values detailed in the function PrintDetail in WormAgent.cpp for every time step of the simulation. The second option is to run several simulations (#define chem_index) and save only the last time point of the simulation to compute the chemotaxis index afterwards (with the chemotaxis_index.ipynb notebook). In this case, also a txt file is generated from the PrintDetail function.
 
The NaCl concentration during cultivation can be changed with the parameter "phenotype(25)" (C_breed). We explored the values C_breed=25,50,100.

### WormAgent.h
WormAgent.h
You can choose between a step of NaCl (#define GRAD_STEP), or the gradient use in the paper (#define GRAD_GAUS). Beware that if you use the step concentration, you need to change the constant Preexposure to a value lower than 5000 for memory reasons.



