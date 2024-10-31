// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (pimple.momentumPredictor())
{
	solve(UEqn == -fvc::grad(p));
    // eqnResidualU = solve(UEqn == -fvc::grad(p)).initialResidual();
    // eqnResidualUx = eqnResidualU[0];
    // eqnResidualUy = eqnResidualU[1];
    // eqnResidualUz = eqnResidualU[2];
}
// else
// {
//     eqnResidualUx = 0;
//     eqnResidualUy = 0;
//     eqnResidualUz = 0;
// }
