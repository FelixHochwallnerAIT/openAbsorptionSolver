    // initialize values for convergence checks

    scalar eqnResidualT = 1; 
    scalar convergenceCriterionT = 0;

    scalar eqnResidualCw = 1; 
    scalar convergenceCriterionCw = 0;

    residualControl.readIfPresent("T", convergenceCriterionT);
    residualControl.readIfPresent("cw", convergenceCriterionCw);