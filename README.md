learning-and-adaptation
=======================

Gaussian Process Optimization based Learning for Trajectory Tracking, Okan Koc, Gajamohan Mohanarajah and Andreas Krause

=======================

TODO List for TGP (Tracking with Gaussian Processes):

- Implement the saturating cost and compare performance with quadratic cost
- Make simple regret experiments to compare the acquisition functions of GPUCB, EI and the new algorithm GP-MI
- Include acquisition functions as a subclass of contextual bandits
- Find an implemented version of PILCO + reference tracking + nominal model (Multi-Task PILCO, P.Englert)
- Learn faster in complicated dynamical mismatches. Things to try:

  1. reward shaping, 
  2. feedback added learning (indirect model learning)
  3. conditioning on estimation data
  4. using options in a hierarchical bandits setting [can mpc be incorporated to this approach with predetermined/flexible horizon?]
  5. parametrize inputs cleverly [Gaussians or time varying linear feedback control structure?]
  6. oracles: phasing as in DMP to get smooth approximating trajectories [could parameters be optimized via RKSH norm of cost differences?]

Remarks:

- Can one smoothen inputs by penalizing input exploration and still achieve no-regret?
- For finite horizon problems, it makes sense to explore progressively towards the end (cautious exploration)
- TGP with robust trajectory generation as a point tracking algorithm to compare with PILCO
