/*---------------------------------------------------------------------------*\
Description
    Calculates the water vapor pressure pw of the LiBr-water solution
    using Patterson's correlation (1984).
\*---------------------------------------------------------------------------*/

// Check if in region solution
if (fluidRegions[i].name().find("solution") == std::string::npos)
{
    Info << "\nError in pwSolution: Region where pwSolution "
         << "should be evaluated is no solution region." 
         << nl << endl;
    exit(1);
}
else
{
    // Patterson's correlation coefficients
    RectangularMatrix<scalar> B(6, 3, 0);
    B[0][0] = -1.313448e-1;  B[0][1] =  9.967944e-1;  B[0][2] =  1.978788e-5;
    B[1][0] =  1.820914e-1;  B[1][1] =  1.778069e-3;  B[1][2] = -1.779481e-5;
    B[2][0] = -5.177356e-2;  B[2][1] = -2.215597e-4;  B[2][2] =  2.002427e-6;
    B[3][0] =  2.827426e-3;  B[3][1] =  5.913618e-6;  B[3][2] = -7.667546e-8;
    B[4][0] = -6.380541e-5;  B[4][1] = -7.308556e-8;  B[4][2] =  1.201525e-9;
    B[5][0] =  4.340498e-7;  B[5][1] =  2.788472e-10; B[5][2] = -6.641716e-12;

	// T unit
	const dimensionedScalar TUnit
	(
		"TUnit",
		dimensionSet(0,0,0,1,0,0,0),
		scalar(1)
	);

    // p unit
	const dimensionedScalar pUnit
	(
		"pUnit",
		dimensionSet(1,-1,-2,0,0,0,0),
		scalar(1.0)
	);

    // Solution parameters
    // Salt concentration
    const volScalarField C1 = 1 - cw; 
    // Mass fraction in percent
    const volScalarField X_percent = C1 * 100; 
    // Convert K to °F
    const volScalarField T_degF = (T/TUnit) * 9.0 / 5.0 - 459.67; 

    // Calculate dew point temperature (T_dp) using Patterson's correlation
    // Initialize T_dp_degF as a uniform dimensonless field with a value of 0
    volScalarField T_dp_degF
    (
        IOobject
        (
            "T_dp_degF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            T_dp_degF += 
                B[i][j] 
                * Foam::pow(X_percent, i) 
                * Foam::pow(T_degF, j);
        }
    }
    // Convert °F to K
    const volScalarField T_dp_K = (T_dp_degF - 32.0) * 5.0 / 9.0 + 273.15; 

    // Calculate vapor pressure (p_w)
    const volScalarField log_p = 7.05 - (1603.54 / T_dp_K) 
                                       - (104095.5 / Foam::pow(T_dp_K, 2));
    // Convert kPa to Pa
    const volScalarField p_Pa = Foam::pow(10.0, log_p) * 1000; 

    // Assign the result to the output variable
    pw = p_Pa * pUnit;
} 
// ************************************************************************* //
