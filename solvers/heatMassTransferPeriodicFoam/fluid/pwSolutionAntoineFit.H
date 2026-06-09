/*---------------------------------------------------------------------------*\
Description
	Calculates the water vapor pressure pw of the solution 
	according to the Antoine Equation log10(pw-D) = A-B/(C+T),
	where the coefficients A, B and C are fit by a third 
	order polynomial in dependency on the solution 
	concentration X.  
	(i.e. A = k0A + k1A * X + k2A * X^2 + k3A * X^3)
	The coefficients for the fit are searched in the 
	thermophysical properties dictionary in the constant
	folder of the solution. They should be located in the
	sub dictionary mixture/vaporPressureFit. They have to 
	have the proper dimensions.
\*---------------------------------------------------------------------------*/

// Check if in region solution
if (fluidRegions[i].name().find("solution") == std::string::npos)
{
	Info<< "\nError in pwSolution: Region where pwSolution "
		<< "should be evaluated is no solution region." 
		<< nl << endl;
	exit(1);
}

else
{
	// Get fitting parameters
	
	// A-fitting parameters
	const dimensionedScalar& k0A
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k0A")
	);
	const dimensionedScalar& k1A
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k1A")
	);
	const dimensionedScalar& k2A
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k2A")
	);
	const dimensionedScalar& k3A
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k3A")
	);

	// B-fitting parameters
	const dimensionedScalar& k0B
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k0B")
	);
	const dimensionedScalar& k1B
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k1B")
	);
	const dimensionedScalar& k2B
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k2B")
	);
	const dimensionedScalar& k3B
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k3B")
	);

	// C-fitting parameters
	const dimensionedScalar& k0C
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k0C")
	);
	const dimensionedScalar& k1C
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k1C")
	);
	const dimensionedScalar& k2C
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k2C")
	);
	const dimensionedScalar& k3C
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("k3C")
	);

	// D-fitting parameter
	const dimensionedScalar& D
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("vaporPressureFit")
		.lookup("D")
	);

	// Get dimensional pressure variable
	const dimensionedScalar p1
	(
		"p1",
		dimensionSet(1,-1,-2,0,0,0,0),
		scalar(1.0)
	);

	// Get dimensional temperature variable
	const dimensionedScalar Tk
	(
		"Tk",
		dimensionSet(0,0,0,1,0,0,0),
		scalar(273.15)
	);

	// pw fit
	// pw = (10^(A(X)-B(X)/(C(X)+T))+D)
	// A(X) = k0A + k1A * X + k2A * X^2 + k3A * X^3
	// B(X) = k0B + k1B * X + k2B * X^2 + k3B * X^3
	// C(X) = k0C + k1C * X + k2C * X^2 + k3C * X^3
	pw = 	(pow(10,
			(k0A+k1A*(1-cw)+k2A*pow((1-cw),2)+k3A*pow((1-cw),3))-
			(k0B+k1B*(1-cw)+k2B*pow((1-cw),2)+k3B*pow((1-cw),3))/
			(k0C+k1C*(1-cw)+k2C*pow((1-cw),2)+k3C*pow((1-cw),3)+
			(T))
			)*p1+D);

}	// End else statement
// ************************************************************************* //
