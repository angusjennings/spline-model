# spline-model
A spline-based approach to smoothly constrain hazard ratios with a view to apply treatment effect waning,
Angus C. Jennings, Mark J. Rutherford, Paul C. Lambert; 
Code used for the paper, 'A spline-based approach to smoothly constrain hazard ratios with a view to apply treatment effect waning', to enforce smooth waning under a flexible parametric model framework. Example and simulation study supporting code.
- ConstrainNR.R; adapted survPen code to allow parameter constraint
- Examples - colon.do; Stata code to reproduce equvalent of example used in the paper (using stpm3)
- Examples - colon.R; R code to reproduce the examples used in the paper (using modified survPen code per ConstrainNR.R)
- Simulation Analyses.R; code used to run the simulation study described (using modified survPen code per ConstrainNR.R)
