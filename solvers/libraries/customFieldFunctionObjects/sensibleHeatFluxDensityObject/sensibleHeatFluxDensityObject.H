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
    Foam::functionObjects::sensibleHeatFluxDensityObject

Description
    Calculates and writes the absorbed mass flux at interface.

    Example of function object specification:
    \verbatim
    sensibleHeatFluxDensityObject1
    {
        type        sensibleHeatFluxDensityObject;
        libs        ("customFieldFunctionObjects.so");
        ...
        patches     (Interface);
    }
    \endverbatim

Usage
    \table
        Property | Description                 | Required    | Default value
        type     | type name: sensibleHeatFluxDensityObject | yes         |
        patches  | list of patches to process  | yes         | 
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    sensibleHeatFluxDensityObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sensibleHeatFluxDensityObject_H
#define functionObjects_sensibleHeatFluxDensityObject_H

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
                       Class sensibleHeatFluxDensityObject Declaration
\*---------------------------------------------------------------------------*/

class sensibleHeatFluxDensityObject
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
        tmp<volScalarField> calcSensibleHeatFluxDensity
		(
			const volScalarField& T
		);

public:

    //- Runtime type information
    TypeName("sensibleHeatFluxDensityObject");


    // Constructors

        //- Construct from Time and dictionary
        sensibleHeatFluxDensityObject
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        sensibleHeatFluxDensityObject(const sensibleHeatFluxDensityObject&) = delete;


    //- Destructor
    virtual ~sensibleHeatFluxDensityObject();


    // Member Functions

        //- Read the sensibleHeatFluxDensityObject data
        virtual bool read(const dictionary&);

        //- Calculate the absorbed mass flux
        virtual bool execute();

        //- Write the absorbed mass flux
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sensibleHeatFluxDensityObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
