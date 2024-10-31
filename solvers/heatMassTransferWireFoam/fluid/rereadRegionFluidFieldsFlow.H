volVectorField U = volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        fluidRegions[i]
    );

phiFluid[i] = surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );