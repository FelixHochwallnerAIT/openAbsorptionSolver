// Calculate phi
volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
	"phiHbyA",
    fvc::flux(HbyA)
    + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU);

// Pimple Non Orthogonal corrector loops
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
    	fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();
    // eqnResidualP = pEqn.solve().initialResidual();

    if (pimple.finalNonOrthogonalIter())
    {
    	// Calculate the conservative fluxes
        phi = phiHbyA - pEqn.flux();
    }

    #include "writeProbes.H"
}



// Solve continuity
#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
