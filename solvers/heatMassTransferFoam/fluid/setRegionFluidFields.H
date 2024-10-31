    // Set fluid mesh
    const fvMesh& mesh = fluidRegions[i];

    // Set fluid field
    surfaceScalarField& phi = phiFluid[i];

    // Set thermophysical properties dict
    IOdictionary& thermophysicalProperties = 
            thermophysicalPropertiesDictFluid[i];

    // Set humidity field
    volScalarField& cw = cwFluid[i];
    dimensionedScalar& cwDiff = cwDiffFluid[i];
    volScalarField& pw = pwFluid[i];

    // Set temperature field
    volScalarField& T = TFluid[i];
    dimensionedScalar& TDiff = TDiffFluid[i];
