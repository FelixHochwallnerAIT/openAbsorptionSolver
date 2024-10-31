/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "sensibleHeatFluxScalarTransport.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sensibleHeatFluxScalarTransport, 0);
    addToRunTimeSelectionTable(functionObject, sensibleHeatFluxScalarTransport, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::sensibleHeatFluxScalarTransport::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Sensible Heat Flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
	writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::sensibleHeatFluxScalarTransport::calcSensibleHeatFlux
(
	const volScalarField& T
)
{

	// Read thermophysicalProperties Dict
	const dictionary& transportProperties =
		lookupObject<IOdictionary>("transportProperties");

	// Extract lambda
	const dimensionedScalar& lambda
	(
		transportProperties
		.lookup("lambda")
	);

    tmp<volScalarField> tsensibleHeatFluxScalarTransport
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar("Q_dot",dimensionSet(1,2,-3,0,0,0,0),scalar(0))
        )
    );

    volScalarField::Boundary& sensibleHeatFluxScalarTransportBf =
        tsensibleHeatFluxScalarTransport.ref().boundaryFieldRef();

	const volScalarField::Boundary& TBf = 
		T.boundaryField();

	const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

    	sensibleHeatFluxScalarTransportBf[patchi] = lambda.value()*TBf[patchi].snGrad()*magSf[patchi]*(-1); 
	}

    return tsensibleHeatFluxScalarTransport;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sensibleHeatFluxScalarTransport::sensibleHeatFluxScalarTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sensibleHeatFluxScalarTransport::~sensibleHeatFluxScalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sensibleHeatFluxScalarTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    return true;
}


bool Foam::functionObjects::sensibleHeatFluxScalarTransport::execute()
{
	const tmp<volScalarField>& T = mesh_.lookupObject<volScalarField>("T");

    word name(type());

    return store(name, calcSensibleHeatFlux(T));
}


bool Foam::functionObjects::sensibleHeatFluxScalarTransport::write()
{
    Log << type() << " " << name() << " write:" << nl;

    //writeLocalObjects::write();

    logFiles::write();

    const volScalarField& sensibleHeatFluxScalarTransport =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& mfp = sensibleHeatFluxScalarTransport.boundaryField()[patchi];

        const scalar minMfp = gMin(mfp);
        const scalar maxMfp = gMax(mfp);
        const scalar integralMfp = gSum(mfp);

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minMfp
                << tab << maxMfp
                << tab << integralMfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minMfp << ", " << maxMfp << ", " << integralMfp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
