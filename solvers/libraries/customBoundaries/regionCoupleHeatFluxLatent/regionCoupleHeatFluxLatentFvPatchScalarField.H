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

Class
    regionCoupleHeatFluxLatentFvPatchScalarField

Description
	Region couple boundary condition for heat flux, including the latent
	heat from the absorbed water. Region coupling at the multi-region 
	interface is accomplished via patchToPatchInterpolation, which is 
	non-conservative for non-conformal boundary patches.

    Example usage:
		interfacePatch
        {
            type                	regionCoupleHeatFluxLatent;
            neighbourRegionName 	neighbourRegion;
            neighbourPatchName  	neighbourPatch;
            neighbourFieldName  	T;
            value               	uniform 0;
        }

SourceFiles
    regionCoupleHeatFluxLatentFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleHeatFluxLatentFvPatchScalarField_H
#define regionCoupleHeatFluxLatentFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "coupleManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupleHeatFluxLatentFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleHeatFluxLatentFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Couple manager object
        coupleManager coupleManager_;


public:

    //- Runtime type information
    TypeName("regionCoupleHeatFluxLatent");


    // Constructors

        //- Construct from patch and internal field
        regionCoupleHeatFluxLatentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleHeatFluxLatentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupleHeatFluxLatentFvPatchScalarField
        //  onto a new patch
        regionCoupleHeatFluxLatentFvPatchScalarField
        (
            const regionCoupleHeatFluxLatentFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new regionCoupleHeatFluxLatentFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupleHeatFluxLatentFvPatchScalarField
        (
            const regionCoupleHeatFluxLatentFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleHeatFluxLatentFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Return the maximum normalized coupled patch residual
	scalar maxResidual() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
