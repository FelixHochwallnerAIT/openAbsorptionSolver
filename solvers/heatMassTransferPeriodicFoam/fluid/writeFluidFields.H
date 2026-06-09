    forAll(fluidRegions, i)
    {
        // Set fluid fields
        surfaceScalarField& phi = phiFluid[i];
        volScalarField& cw = cwFluid[i];
        volScalarField& pw = pwFluid[i];
        volScalarField& T = TFluid[i];

        // Write fluid fields
        phi.write();
        cw.write();
        pw.write();
        T.write();
    }