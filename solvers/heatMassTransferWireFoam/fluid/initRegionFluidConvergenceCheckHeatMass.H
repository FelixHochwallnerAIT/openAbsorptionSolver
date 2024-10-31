    // initialize values for convergence checks

    scalar eqnResidualT = 1; 
    scalar convergenceCriterionT = 0;

    scalar eqnResidualCw = 1; 
    scalar convergenceCriterionCw = 0;

    heatMassDict.subDict("residualControl").readIfPresent("T", convergenceCriterionT);
    heatMassDict.subDict("residualControl").readIfPresent("cw", convergenceCriterionCw);
