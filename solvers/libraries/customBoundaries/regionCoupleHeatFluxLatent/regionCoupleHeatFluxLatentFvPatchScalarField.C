/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "regionCoupleHeatFluxLatentFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "regionCoupleTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleHeatFluxLatentFvPatchScalarField::
regionCoupleHeatFluxLatentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p)
{}


Foam::regionCoupleHeatFluxLatentFvPatchScalarField::
regionCoupleHeatFluxLatentFvPatchScalarField
(
    const regionCoupleHeatFluxLatentFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_)
{}


Foam::regionCoupleHeatFluxLatentFvPatchScalarField::
regionCoupleHeatFluxLatentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p, dict)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::regionCoupleHeatFluxLatentFvPatchScalarField::
regionCoupleHeatFluxLatentFvPatchScalarField
(
    const regionCoupleHeatFluxLatentFvPatchScalarField& whftcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(whftcsf, iF),
    coupleManager_(whftcsf.coupleManager_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleHeatFluxLatentFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
	
	// Get neighbour field
	const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

	// Get flux of neighbour field
	const scalarField qNbr = refCast<const regionCoupleTemperatureFvPatchScalarField>(Tnbr).flux();

	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");

	// Read heat conductivity lambda
	const dimensionedScalar& lambda
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("lambda")
	);

	// Read h condensation
	const dimensionedScalar& deltaH
	(
		coupleManager_.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("humidity")
		.lookup("deltaH")
	);

	// Get cw Field of neighbour region
	const fvPatchScalarField& cwA
	(
		coupleManager_
		.neighbourPatch()
		.lookupPatchField
		<GeometricField<scalar, fvPatchField, volMesh>, scalar>
		("cw")
	);

	// Get diffusion coefficient
	const dimensionedScalar& cwDiffA
	(
		coupleManager_
		.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("transport")
		.lookup("cwDiff")
	);
	
	// Get density of neighbour region
	const dimensionedScalar& rhoA
	(
		coupleManager_
		.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("equationOfState")
		.lookup("rho")
	);

	// Calculate mass flux
	const scalarField mDot = rhoA.value()*cwDiffA.value()/(1-cwA)*cwA.snGrad();
	
	// Calculate latent heat flux
	const scalarField qDotLatent = deltaH.value()*mDot;

	// Enforce flux matching
	gradient() = ((qNbr+qDotLatent)/lambda.value())*(-1);

	fixedGradientFvPatchScalarField::updateCoeffs();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleHeatFluxLatentFvPatchScalarField::maxResidual() const
{
	// Get neighbour field
	const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

	// Get flux of neighbour field
	const scalarField qNbr = refCast<const regionCoupleTemperatureFvPatchScalarField>(Tnbr).flux();

	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");

	// Read heat conductivity lambda
	const dimensionedScalar& lambda
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("lambda")
	);

	// Read h condensation
	const dimensionedScalar& deltaH
	(
		coupleManager_.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("humidity")
		.lookup("deltaH")
	);

	// Get cw Field
	const fvPatchScalarField& cwA
	(
		coupleManager_
		.neighbourPatch()
		.lookupPatchField
		<GeometricField<scalar, fvPatchField, volMesh>, scalar>
		("cw")
	);
	// Get density and diffusion
	const dimensionedScalar& cwDiffA
	(
		coupleManager_
		.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("transport")
		.lookup("cwDiff")
	);
	const dimensionedScalar& rhoA
	(
		coupleManager_
		.neighbourRegion()
		.lookupObject<IOdictionary>("thermophysicalProperties")
		.subDict("mixture")
		.subDict("equationOfState")
		.lookup("rho")
	);

	// Calculate mass flux
	const scalarField mDot = rhoA.value()*cwDiffA.value()/(1-cwA)*cwA.snGrad();
	
	// Calculate latent heat flux
	const scalarField qDotLatent = deltaH.value()*mDot;

    // Calculate the maximum normalized residual
    const fvPatchScalarField& Town = *this;
    const scalarField& qOwn = lambda.value()*Town.snGrad();
    const scalar& residual =
	gMax
	(
		mag(mag(qOwn) - mag(qNbr) - mag(qDotLatent))/
		max(min(gMax(mag(qOwn)),gMax(mag(qNbr)+mag(qDotLatent))), SMALL)
	);

    return residual;
}


//- Write
void Foam::regionCoupleHeatFluxLatentFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleHeatFluxLatentFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
