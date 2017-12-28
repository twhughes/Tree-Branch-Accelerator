# DLA-Laser-Coupling-Simulation-Software

## Contents
This software package is used to simulate an optical power delivery system for dielectric laser accelerators.  The structure consists of a fractal waveguide network as described in the accompanying paper.  This code takes as input all of the assumed parameters for the system (defaults assigned automatically) and computes the energy gain and acceleration gradient from the resulting structure.  First, all of the input fields necessary to cause:
- damage at the accelerator structures
- damage at the input facet
- self-phase modulation and pulse degradation in the waveguides
- self-focusing effects in the waveguides
are computed.  The minimum constraint is then set as the limiting input electric field.  Perfectly phased output ports are assumed in the simulations.

## Paper information
Accompanying a paper on on-chip laser coupling for dielectric laser accelerator structures.  Right now it is an arXiv preprint.

[ArXiv preprint](https://arxiv.org/abs/1709.04441 "ArXiv preprint")

## How to run simulation
The main codebase is in "simulation.m"  This defines the simulation class, which holds both the properties of parameters and methods corresponding to solving the laser coupling system.

Simple examples are in the appropriately named ```simple_examples.m```, which will help new users get to know how the code works and solve some simple examples.

For documentation, type:

        doc simulation

in the matlab command line once the directory is added to path

There are auxilary functions in the ```components``` directory, which are used to simulate the DLA structures among other things.

## Citing our work

If you would like to use this code, please cite the paper:


        @article{hughes2017chip,
          title={On-Chip Laser Power Delivery System for Dielectric Laser Accelerators},
          author={Hughes, Tyler W and Tan, Si and Zhao, Zhexin and Sapra, Neil V and Lee, Yun Jo and Leedle, Kenneth J and Deng, Huiyang and Miao, Yu and Black, Dylan S and Qi, Minghao and others},
          journal={arXiv preprint arXiv:1709.04441},
          year={2017}
        }

