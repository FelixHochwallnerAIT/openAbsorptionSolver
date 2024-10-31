/*---------------------------------------------------------------------------*\
Description
	Calculates the water vapor pressure pw of the LiBr solution.
\*---------------------------------------------------------------------------*/

// Check if in region solution
if (fluidRegions[i].name() != "solution")
{
	Info<< "\nError in pwSolution: Region where pwSolution "
		<< "should be evaluated is no solution region." 
		<< nl << endl;
	exit(1);
}

else
{
	// Constants
	const scalar ups=2;
	// molar mass of LiBr
	const scalar M1 = 86.65; // g/mol
	// molar mass of water
	const scalar M2 = 18.015; // g/mol
	// T0 = 275.15 K
	const dimensionedScalar T0
	(
		"T0",
		dimensionSet(0,0,0,1,0,0,0),
		scalar(273.15)
	);
	// Tunit = 1 K
	const dimensionedScalar Tunit
	(
		"Tunit",
		dimensionSet(0,0,0,1,0,0,0),
		scalar(1)
	);

	// a-vector
	RectangularMatrix<scalar> a(1,3,0);
	a[0][0] =  1.1375e1;
	a[0][1] = -3.8590e3;
	a[0][2] =  5.1319e5;

	// b-vector
	RectangularMatrix<scalar> b(1,3,0);
	b[0][0] =  8.6010e-1;
	b[0][1] = -1.9575e2;
	b[0][2] =  2.3136e4;

	// A-matrix
	RectangularMatrix<scalar> A(6,3,0);
	A[0][0] = -2.1963155e1;
	A[0][1] =  4.9372316e3;
	A[0][2] = -6.5548406e5;
	A[1][0] = -3.8104752e3;
	A[1][1] =  2.6115345e6;
	A[1][2] = -3.6699691e8;
	A[2][0] =  1.2280854e5;
	A[2][1] = -7.7187923e7;
	A[2][2] =  1.0398560e10;
	A[3][0] = -1.4716737e6;
	A[3][1] =  9.1952848e8;
	A[3][2] = -1.1894502e11;
	A[4][0] =  7.7658213e6;
	A[4][1] = -4.9375666e9;
	A[4][2] =  6.3175547e11;
	A[5][0] = -1.5118922e7;
	A[5][1] =  9.8399744e9;
	A[5][2] = -1.2737898e12;

	// Salt concentration
	const volScalarField C1 = 1-cw;

	// Mean molar mass
	const volScalarField Mm = 1/(C1/M1+(1-C1)/M2);

	const volScalarField x  = C1*Mm/M1;
	const volScalarField m  = x/((1-x)*M2);

	const volScalarField Tul = T/Tunit;

	const volScalarField alpha = a[0][0] + a[0][1]/Tul + a[0][2]/pow(Tul,2);
	const volScalarField beta  = b[0][0] + b[0][1]/Tul + b[0][2]/pow(Tul,2);

	const volScalarField phi
	(
		1 +
		(A[0][0] + A[0][1]/Tul + A[0][2]/pow(Tul,2))*pow(m,0.5)+
		(A[1][0] + A[1][1]/Tul + A[1][2]/pow(Tul,2))*pow(m,1.0)+
		(A[2][0] + A[2][1]/Tul + A[2][2]/pow(Tul,2))*pow(m,1.5)+
		(A[3][0] + A[3][1]/Tul + A[3][2]/pow(Tul,2))*pow(m,2.0)+
		(A[4][0] + A[4][1]/Tul + A[4][2]/pow(Tul,2))*pow(m,2.5)+
		(A[5][0] + A[5][1]/Tul + A[5][2]/pow(Tul,2))*pow(m,3.0)
	);
	
	const volScalarField Tt = (T-T0)/Tunit;

	const volScalarField pStar = 0.611*Foam::exp(17.62*Tt/(243.12+Tt));

	const volScalarField theta = Foam::cosh(alpha-beta*Foam::log(pStar))*Foam::exp(phi*ups*m*M2*beta);

	const volScalarField pSol = Foam::exp(1/beta*(alpha-Foam::log(theta+pow(pow(theta,2)-1,0.5))));

	// Get dimensional pressure variable
	const dimensionedScalar p1
	(
		"p1",
		dimensionSet(1,-1,-2,0,0,0,0),
		scalar(1.0)
	);

	// Get dimensional partial vapor pressure
	pw = pSol*p1*1000;

}	// End else statement
// ************************************************************************* //
