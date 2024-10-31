    // Initialise fluid field pointer lists

    // Flow field
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
        
    // Thermophysical properties dict
    PtrList<IOdictionary> thermophysicalPropertiesDictFluid(fluidRegions.size());

    // Humidity field
    PtrList<dimensionedScalar> cwDiffFluid(fluidRegions.size());
    PtrList<volScalarField> cwFluid(fluidRegions.size());
    PtrList<volScalarField> pwFluid(fluidRegions.size());

    // Temperature field
    PtrList<dimensionedScalar> TDiffFluid(fluidRegions.size());
    PtrList<volScalarField> TFluid(fluidRegions.size());


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {

        Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

        Info<< "    Reading U\n" << endl;
        volVectorField U
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

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
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
            )
        );


        // Thermophysical properties dict
        Info<< "    Adding to thermophysicalPropertiesDictFluid\n" << endl;
        // Thermophysical properties of fluid as defined in dict
        thermophysicalPropertiesDictFluid.set 
        (
                i,
                new IOdictionary
                (
                    IOobject
                    (
                        "thermophysicalProperties",
                        runTime.constant(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
        );

                
        // Humidity field
        Info<< "    Adding to cwDiffFluid\n" << endl;
        cwDiffFluid.set // Diffusion of water in air
        (
            i,
                new dimensionedScalar
                (
                    thermophysicalPropertiesDictFluid[i]
                    .subDict("mixture")
                    .subDict("transport")
                    .lookup("cwDiff")
                )
        );

        Info<< "    Adding to cwFluid\n" << endl;
        cwFluid.set // Water mass fraction
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "cw",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to pwFluid\n" << endl;
        pwFluid.set // Water vapor pressure
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "pw",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
            fluidRegions[i]
            )
        );


        // Temperature field
        Info<< "    Adding to lambdaFluid\n" << endl;
        TDiffFluid.set // Diffusion of temperature
        (
            i,
            new dimensionedScalar
            (
                thermophysicalPropertiesDictFluid[i]
                .subDict("mixture")
                .subDict("transport")
                .lookup("TDiff")
            )
        );
        
        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set // Temperature
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
            fluidRegions[i]
            )
        );



    }


