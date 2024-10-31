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

#include "regionCoupleVaporPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleVaporPressureFvPatchScalarField::
regionCoupleVaporPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p),
    relax_(0.5)
{}


Foam::regionCoupleVaporPressureFvPatchScalarField::
regionCoupleVaporPressureFvPatchScalarField
(
    const regionCoupleVaporPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    relax_(ptf.relax_)
{}


Foam::regionCoupleVaporPressureFvPatchScalarField::
regionCoupleVaporPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p, dict),
    relax_(readScalar(dict.lookup("relaxationFactor")))
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


Foam::regionCoupleVaporPressureFvPatchScalarField::
regionCoupleVaporPressureFvPatchScalarField
(
    const regionCoupleVaporPressureFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    relax_(wtcsf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleVaporPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get neighbour field
	const fvPatchField<scalar>& pwNbr = coupleManager_.neighbourPatchField<scalar>();

	// Get p0 and beta from thermophysicalProperties
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");

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

	const fvPatchField<scalar>& pwOwn = patch().lookupPatchField<volScalarField, scalar>("pw");

	// Relax vapor pressure equilibrium
	const Field<scalar> pwRelax = pwOwn + relax_*(pwNbr - pwOwn);
	
	// Enforce vapor pressure equilibrium
	operator==(beta.value()*pwRelax/(p0.value()-((1-beta.value())*pwRelax)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


//- Return the patch flux
Foam::tmp<Foam::scalarField>
Foam::regionCoupleVaporPressureFvPatchScalarField::flux() const
{
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

	// Scalar field cw
	const fvPatchScalarField& cw = *this;

	// Eckert Schneider dotM/A = rho*D/(1-cw)*dcw/dy
	return rho.value()*cwDiff.value()/(1-cw)*cw.snGrad();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleVaporPressureFvPatchScalarField::maxResidual() const
{
	// Get neighbour field
	const fvPatchField<scalar>& pwNbr = coupleManager_.neighbourPatchField<scalar>();

    // Calculate the maximum normalized residual
	const fvPatchField<scalar>& pwOwn = patch().lookupPatchField<volScalarField, scalar>("pw");
    const scalar& residual =
        gMax
        (
            mag(pwOwn - pwNbr)/
            max(min(gMax(pwOwn),gMax(pwNbr)), SMALL)
        );
        
    return residual; 
}


//- Write
void Foam::regionCoupleVaporPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    os.writeKeyword("relaxationFactor") << relax_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleVaporPressureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
