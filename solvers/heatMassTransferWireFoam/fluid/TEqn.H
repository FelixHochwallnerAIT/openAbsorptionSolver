// Convection-Diffusion Equation
fvScalarMatrix TEqn
(
    fvm::ddt(T)
  + fvm::div(phi, T)
  - fvm::laplacian(TDiff, T)
);

TEqn.relax();
eqnResidualT = TEqn.solve().initialResidual();

Info<< "Min/max T:" << min(T).value() << ' '
	<< max(T).value() << endl;
