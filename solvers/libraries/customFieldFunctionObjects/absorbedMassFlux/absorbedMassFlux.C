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

#include "absorbedMassFlux.H"
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
    defineTypeNameAndDebug(absorbedMassFlux, 0);
    addToRunTimeSelectionTable(functionObject, absorbedMassFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::absorbedMassFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Absorbed Mass Flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
	writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::absorbedMassFlux::calcMassFlux
(
	const volScalarField& cw
)
{

	// Read thermophysicalProperties Dict
	const dictionary& thermophysicalProperties =
		lookupObject<IOdictionary>("thermophysicalProperties");

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

    tmp<volScalarField> tabsorbedMassFlux
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar("m_abs",dimensionSet(1,0,-1,0,0,0,0),scalar(0))
        )
    );

    volScalarField::Boundary& absorbedMassFluxBf =
        tabsorbedMassFlux.ref().boundaryFieldRef();

	const volScalarField::Boundary& cwBf = 
		cw.boundaryField();

	const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

    	absorbedMassFluxBf[patchi] = rho.value()*cwDiff.value()/(1-cwBf[patchi])
									*cwBf[patchi].snGrad()*magSf[patchi]*(-1); 
	}

    return tabsorbedMassFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::absorbedMassFlux::absorbedMassFlux
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

Foam::functionObjects::absorbedMassFlux::~absorbedMassFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::absorbedMassFlux::read(const dictionary& dict)
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


bool Foam::functionObjects::absorbedMassFlux::execute()
{
	const tmp<volScalarField>& cw = mesh_.lookupObject<volScalarField>("cw");

    word name(type());

    return store(name, calcMassFlux(cw));
}


bool Foam::functionObjects::absorbedMassFlux::write()
{
    Log << type() << " " << name() << " write:" << nl;

    //writeLocalObjects::write();

    logFiles::write();

    const volScalarField& absorbedMassFlux =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& mfp = absorbedMassFlux.boundaryField()[patchi];

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
