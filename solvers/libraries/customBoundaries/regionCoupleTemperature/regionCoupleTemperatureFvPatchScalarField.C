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

#include "regionCoupleTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p),
    relax_(0.5)
{}


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const regionCoupleTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    relax_(ptf.relax_)
{}


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
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


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const regionCoupleTemperatureFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    relax_(wtcsf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get neighbour field
	const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

	// Enforce temperature boundary condition
	operator==(*this + relax_*(Tnbr - *this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


//- Return the patch flux
Foam::tmp<Foam::scalarField>
Foam::regionCoupleTemperatureFvPatchScalarField::flux() const
{
	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		db().lookupObject<IOdictionary>("thermophysicalProperties");
	
	// Get heat conduction lambda
	const dimensionedScalar& lambda
	(
		thermophysicalProperties
		.subDict("mixture")
		.subDict("transport")
		.lookup("lambda")
	);

	// Scalar field T
	const fvPatchScalarField& T = *this;

	// Thermal conduction: heat flux density q = lambda*dT/dy
	return lambda.value()*T.snGrad();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleTemperatureFvPatchScalarField::maxResidual() const
{
    // Get neighbour field
    const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

    // Calculate the maximum normalized residual
    const scalarField& Town = *this;
    const scalar& residual =
        gMax
        (
            mag(Town - Tnbr)/
            max(min(gMax(Town),gMax(Tnbr)), SMALL)
        );
        
    return residual; 
}


//- Write
void Foam::regionCoupleTemperatureFvPatchScalarField::write
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
    regionCoupleTemperatureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
