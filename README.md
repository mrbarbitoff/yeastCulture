# Yeasst \[PSI+\] prion propagation model.
This repository contains the recreated \[PSI+\] model from *Derdowski et al.* 

Due to certain inconsistewncies with the original findings reoirted by Derdowski et al., the model was implemented using several languages:
* MATLAB (original implementation, 2015);
* R (2015)
* Python (2015)
* **In 2023, the model was re-implemented using Kotlin**

The current Kotlin implementation includes a `YeastCell` class, for which several custom extension functions were implemented:
* `YeastCell.divide()` - executes cell division in accordance with the rules from the original study, return as `YeastCell1 object (daughter cell);
* `YeastCell.nextReaction()` - simulates the next reaction (conversion, fragmentation, or protein synthesis), and its respective time interval from the exponential distribution;
* `YeastCell.updateProteinLevel()` - increases the proten level up until the next division in cells that lost the prion.
* `YeastCell.runSimulation` - runs the simulation for the specified time pariod `simTime`, performing reactions and cell divisions. Returns a list of `YeastCell` objects.

The `main` function executes a simulation for a range of $\beta$ and $\gamma$ values. The function uses a `.parmap()` generic extension taken from the `kutils` package.