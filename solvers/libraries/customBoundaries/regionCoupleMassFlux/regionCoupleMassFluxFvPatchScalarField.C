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

#include "regionCoupleMassFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "regionCoupleVaporPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleMassFluxFvPatchScalarField::
regionCoupleMassFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p)
{}


Foam::regionCoupleMassFluxFvPatchScalarField::
regionCoupleMassFluxFvPatchScalarField
(
    const regionCoupleMassFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_)
{}


Foam::regionCoupleMassFluxFvPatchScalarField::
regionCoupleMassFluxFvPatchScalarField
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


Foam::regionCoupleMassFluxFvPatchScalarField::
regionCoupleMassFluxFvPatchScalarField
(
    const regionCoupleMassFluxFvPatchScalarField& whftcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(whftcsf, iF),
    coupleManager_(whftcsf.coupleManager_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleMassFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get neighbour field
	const fvPatchScalarField& cwNbr = coupleManager_.neighbourPatchField<scalar>();

	const scalarField mDotNbr
	(
		refCast<const regionCoupleVaporPressureFvPatchScalarField>
		(cwNbr).flux()
	);

	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");

	// Extract density and diffusion
	const dimensionedScalar& cwDiff
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("cwDiff")
	);
	const dimensionedScalar& rho
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("equationOfState")
		.lookup("rho")
	);

	// Initalize scalar field cw 
	// from previous time step (first guess)
	//fvPatchScalarField cw = *this;
	scalarField cw = *this;

	// Initalize scalar field cwOld
	// for residual calculation
	// scalarField cwOld = cw;

	// Initialize gradient
	scalarField grad
	(
		(mDotNbr*(1-cw)/(rho.value()*cwDiff.value()))*(-1)
	);
	// scalarField grad
	// (
	// 	(mDotNbr*(1-this->patchInternalField())/
    //     (rho.value()*cwDiff.value()))*(-1)
	// );
	// Print max and min of gradient
	Info<< "Min/max gradient: " << min(grad) << " / "
		<< max(grad) << endl;

	gradient() = grad;
	// Only absorption possible
	//gradient() = max(scalarField(0*cw), grad);

	// Enforce flux matching
	fixedGradientFvPatchScalarField::updateCoeffs();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleMassFluxFvPatchScalarField::maxResidual() const
{
    // Get neighbour field
    const fvPatchField<scalar>& cwNbr = coupleManager_.neighbourPatchField<scalar>();

    scalarField mDotNbr =
        (
            refCast<const regionCoupleVaporPressureFvPatchScalarField>
            (cwNbr).flux()
        );

	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");

	// Extract density and diffusion
	const dimensionedScalar& cwDiff
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("cwDiff")
	);
	const dimensionedScalar& rho
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("equationOfState")
		.lookup("rho")
	);

    // Calculate the maximum normalized residual
    const fvPatchScalarField& cwOwn = *this;
    const scalarField& mDotOwn = rho.value()*cwDiff.value()/(1-cwOwn)*cwOwn.snGrad();
    const scalar& residual =
	gMax
	(
		mag(mag(mDotOwn) - mag(mDotNbr))/
		max(min(gMax(mag(mDotOwn)),gMax(mag(mDotNbr))), SMALL)
	);

    return residual;
}


//- Write
void Foam::regionCoupleMassFluxFvPatchScalarField::write
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
    regionCoupleMassFluxFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
