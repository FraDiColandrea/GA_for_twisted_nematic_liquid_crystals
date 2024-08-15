# GA_for_twisted_nematic_liquid_crystals
A genetic algorithm for the response of twisted nematic liquid crystals to an applied field

By Alicia Sit and Francesco Di Colandrea of the University of Ottawa, 2024.
See http://arxiv.org/abs/2403.00948 for more details. For the published manuscript: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.109.054705

# FUNCTIONS
GA.m : Genetic algorithm to numerically calculate integration parameter (beta) and maximum tilt angle (thetaM) for a given maximum twist angle and applied voltage.
 
GAtilt.m: Genetic algorithm to numerically calculate the liquid crystal tilt distribution for a given beta and thetaM.

twist.m: Integral to analytically calculate the tilt angle at position z given a tilt angle at position z.

RunGA.m: Main file to run GA.m. Can set the range of voltages to use, as well as the max twist angle.

RunTiltTwist.m: Main file to run to obtain the tilt and twist distributions.

==================
# Notes:
-Currently, the material parameters for 6CHBT nematic liquid crystals are used.
