    PtrList<fvMesh> fluidRegions(rp["fluid"].size());
    PtrList<word> solveRegionFluidFlow(rp["solveFluidFlow"].size());

    forAll(rp["fluid"], i)
    {
        Info<< "Create fluid mesh for region " << rp["fluid"][i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp["fluid"][i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
        
        solveRegionFluidFlow.set
        (
                i,
                new word(rp["solveFluidFlow"][i])
        );
    }
