# DLA-Laser-Coupling-Simulation-Software

![alt text](https://github.com/twhughes/DLA-Laser-Coupling-Simulation-Software/blob/master/images/structure.png 'test')

## Contents
This software package is used to simulate an optical power delivery system for dielectric laser accelerators.  The structure consists of a fractal waveguide network as described in the accompanying paper.  This code takes as input all of the assumed parameters for the system (defaults assigned automatically) and computes the energy gain and acceleration gradient from the resulting structure with the following procedure:
1. Find the largest input power (given the parameters) to avoid all of the following
- damage at the accelerator structures
- damage at the input facet
- self-phase modulation and pulse degradation in the waveguides
- self-focusing effects in the waveguides
2. Propagate this pulse through the waveguides and splits, incorporating loss.
3. Compute the energy gain in the DLA assuming perfect phasing of output ports.

A plot of several combinations of pulse duration and Q factor is shown here

![alt text](https://github.com/twhughes/DLA-Laser-Coupling-Simulation-Software/blob/master/images/results.png 'test')

## Paper information
This is the code accompanying a paper on on-chip laser coupling for dielectric laser accelerator structures.

[On-Chip Laser Power Delivery System for Dielectric Laser Accelerators](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.9.054017)

## How to run simulation
The main codebase is in "simulation.m"  This defines the simulation class, which holds both the properties of parameters and methods corresponding to solving the laser coupling system.

Simple examples are in the appropriately named ```simple_examples.m```, which will help new users get to know how the code works and solve some simple examples.

For documentation, type:

        doc simulation

in the matlab command line once the directory is added to path

There are auxilary functions in the ```components``` directory, which are used to simulate the DLA structures among other things.

## Citing our work

If you would like to use this code, please cite the paper:


    @article{hughes2018chip,
      title={On-Chip Laser-Power Delivery System for Dielectric Laser Accelerators},
      author={Hughes, Tyler W and Tan, Si and Zhao, Zhexin and Sapra, Neil V and Leedle, Kenneth J and Deng, Huiyang and Miao, Yu and Black, Dylan S and Solgaard, Olav and Harris, James S and others},
      journal={Physical Review Applied},
      volume={9},
      number={5},
      pages={054017},
      year={2018},
      publisher={APS}
    }

