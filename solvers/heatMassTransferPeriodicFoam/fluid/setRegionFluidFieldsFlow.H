    // Set fluid mesh
    const fvMesh& mesh = fluidRegions[i];

    // Set flow field
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& p = pFluid[i];
    dimensionedScalar& nu = nuFluid[i];

    // Set pRef
    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    // Pimple controls
    pimpleNoLoopControl& pimple = pimples.pimple(i);
    scalar cumulativeContErr = cumulativeContErrs[i];
