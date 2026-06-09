volVectorField U = UPeriodic[i_periodic];

UFluid[i] = U;

phiFluid[i] = surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );