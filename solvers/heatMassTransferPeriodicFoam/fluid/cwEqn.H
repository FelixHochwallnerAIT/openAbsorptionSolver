// Convection-Diffusion Equation
fvScalarMatrix cwEqn
(
	fvm::ddt(cw)
  + fvm::div(phi, cw)
  - fvm::laplacian(cwDiff, cw)
);

cwEqn.relax();
eqnResidualCw = cwEqn.solve().initialResidual();
	
// Solving for pw
if (fluidRegions[i].name() == "air")
{
	const dimensionedScalar& p0
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("p0")
	);

	const dimensionedScalar& beta
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("beta")
	);

	pw = cw/(beta+(1-beta)*cw)*p0;
}
else if (fluidRegions[i].name().find("solution") != std::string::npos)
{
	Info<< "Min/max X: " << min(1-cw).value()
		<< " / " << max(1-cw).value() << endl;

	const int& flagMedium
	(readScalar(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("flagMedium")
	));

    if (flagMedium == 0)
    {
        // LiBr
        #include "pwSolutionLiBrPatterson1984.H"
    }
    else if (flagMedium == 1)
    {
        // Using custom Antoine fit
        // pw = (10^(A(X)-B(X)/(C(X)+T))+D)
        // A(X) = k0A + k1A * X + k2A * X^2 + k3A * X^3
        // B(X) = k0B + k1B * X + k2B * X^2 + k3B * X^3
        // C(X) = k0C + k1C * X + k2C * X^2 + k3C * X^3
        #include "pwSolutionAntoineFit.H"
    }
	else
	{
		// No media selected
		const dimensionedScalar& p0
		(
			thermophysicalProperties
			.subDict("mixture")
			.subDict("transport")
			.lookup("p0")
		);
			pw = p0;
	}
}
else
{
	const dimensionedScalar& p0
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("p0")
	);
		pw = p0;
}
Info<< "Min/max cw: " << min(cw).value() << " / "
	<< max(cw).value() << endl;
