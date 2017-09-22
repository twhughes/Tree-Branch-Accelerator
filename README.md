# DLA-Laser-Coupling-Simulation-Software

Accompanying a paper on on-chip laser coupling for dielectric laser accelerator structures.  Right now it is an arXiv preprint.

        https://arxiv.org/abs/1709.04441

The main codebase is in simulation.m  This defines the simulation class, which has properties of parameters and methods corresponding to solving different quantities.

Simple examples are in simple_examples.m to get to know how the code works and solve some simple examples.

type 


        doc simulation

in the matlab command line once added to path for documentation

If you would like to use this code, please cite the paper:


        @article{hughes2017chip,
          title={On-Chip Laser Power Delivery System for Dielectric Laser Accelerators},
          author={Hughes, Tyler W and Tan, Si and Zhao, Zhexin and Sapra, Neil V and Lee, Yun Jo and Leedle, Kenneth J and Deng, Huiyang and Miao, Yu and Black, Dylan S and Qi, Minghao and others},
          journal={arXiv preprint arXiv:1709.04441},
          year={2017}
        }