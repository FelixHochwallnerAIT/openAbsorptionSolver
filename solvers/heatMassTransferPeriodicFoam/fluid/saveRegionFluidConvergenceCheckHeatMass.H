// check convergence
// Info<< "Equation residual T: " << eqnResidualT
//     << "  convergence criterion: " << convergenceCriterionT
//     << endl;

// Info<< "Equation residual cw: " << eqnResidualCw
//     << "  convergence criterion: " << convergenceCriterionCw
//     << endl;

if ((eqnResidualT < convergenceCriterionT) && (eqnResidualCw < convergenceCriterionCw))
{
    Info<< "Reached convergence criterion in region " 
        << fluidRegions[i].name() << "!" << endl;
    // runTime.writeAndEnd();
    // Info<< "latestTime = " << runTime.timeName() << endl;
    convergedList[i] = true;
}
else
{
    convergedList[i] = false;
}