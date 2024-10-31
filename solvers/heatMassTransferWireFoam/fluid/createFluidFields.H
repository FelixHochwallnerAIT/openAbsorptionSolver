/// Initialise fluid field pointer lists

// Flow problem
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<dimensionedScalar> nuFluid(fluidRegions.size());

// pRef
List<label> pRefCellFluid(fluidRegions.size(),0);
List<scalar> pRefValueFluid(fluidRegions.size(),0.0);

// Heat and mass transfer problem
// Humidity field
PtrList<dimensionedScalar> cwDiffFluid(fluidRegions.size());
PtrList<volScalarField> cwFluid(fluidRegions.size());
PtrList<volScalarField> pwFluid(fluidRegions.size());
// Temperature field
PtrList<dimensionedScalar> TDiffFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());

// Thermophysical properties dict
PtrList<IOdictionary> thermophysicalPropertiesDictFluid(fluidRegions.size());


// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

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
		
	// Flow field
    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

   	Info<< "    Adding to nuFluid\n" << endl;
   	nuFluid.set // Diffusion of momentum in fluid
   	(
       	i,
		new dimensionedScalar
		(
			thermophysicalPropertiesDictFluid[i]
			.subDict("mixture")
			.subDict("transport")
			.lookup("nu")
		)
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
            linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to pFluid\n" << endl;
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i]
        )
    );

	// pRef
    setRefCell
    (
        pFluid[i],
        fluidRegions[i].solutionDict().subDict("PIMPLE"),
        pRefCellFluid[i],
        pRefValueFluid[i]
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
    Info<< "    Adding to TDiffFluid\n" << endl;
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
