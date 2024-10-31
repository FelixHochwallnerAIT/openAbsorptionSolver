/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::absorbedMassFluxDensityObject

Description
    Calculates and writes the absorbed mass flux at interface.

    Example of function object specification:
    \verbatim
    absorbedMassFluxDensity1
    {
        type        absorbedMassFluxDensityObject;
        libs        ("customFieldFunctionObjects.so");
        ...
        patches     (Interface);
    }
    \endverbatim

Usage
    \table
        Property | Description                 			| Required    | Default value
        type     | type name: absorbedMassFluxDensity 	| yes         |
        patches  | list of patches to process  			| yes         | 
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    absorbedMassFluxDensity.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_absorbedMassFluxDensityObject_H
#define functionObjects_absorbedMassFluxDensityObject_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class absorbedMassFluxDensity Declaration
\*---------------------------------------------------------------------------*/

class absorbedMassFluxDensityObject
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the mass flux
        tmp<volScalarField> calcMassFlux
		(
			const volScalarField& cw
		);

public:

    //- Runtime type information
    TypeName("absorbedMassFluxDensityObject");


    // Constructors

        //- Construct from Time and dictionary
        absorbedMassFluxDensityObject
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        absorbedMassFluxDensityObject(const absorbedMassFluxDensityObject&) = delete;


    //- Destructor
    virtual ~absorbedMassFluxDensityObject();


    // Member Functions

        //- Read the absorbedMassFluxDensityObject data
        virtual bool read(const dictionary&);

        //- Calculate the absorbed mass flux
        virtual bool execute();

        //- Write the absorbed mass flux
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const absorbedMassFluxDensityObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
