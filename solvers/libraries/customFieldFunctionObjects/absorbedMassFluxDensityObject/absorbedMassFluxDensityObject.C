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

#include "absorbedMassFluxDensityObject.H"
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
    defineTypeNameAndDebug(absorbedMassFluxDensityObject, 0);
    addToRunTimeSelectionTable(functionObject, absorbedMassFluxDensityObject, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::absorbedMassFluxDensityObject::writeFileHeader(const label i)
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
Foam::functionObjects::absorbedMassFluxDensityObject::calcMassFlux
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

    tmp<volScalarField> tabsorbedMassFluxDensity
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar("m_abs",dimensionSet(1,0,-1,0,0,0,0),scalar(0))
        )
    );

    volScalarField::Boundary& absorbedMassFluxDensityBf =
        tabsorbedMassFluxDensity.ref().boundaryFieldRef();

	const volScalarField::Boundary& cwBf = 
		cw.boundaryField();
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

    	absorbedMassFluxDensityBf[patchi] = rho.value()*cwDiff.value()/(1-cwBf[patchi])
									*cwBf[patchi].snGrad()*(-1); 
	}

    return tabsorbedMassFluxDensity;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::absorbedMassFluxDensityObject::absorbedMassFluxDensityObject
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

Foam::functionObjects::absorbedMassFluxDensityObject::~absorbedMassFluxDensityObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::absorbedMassFluxDensityObject::read(const dictionary& dict)
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


bool Foam::functionObjects::absorbedMassFluxDensityObject::execute()
{

/*  
	tmp<volScalarField> cw;
    if (mesh_.foundObject<volScalarField>("cw"))
    {
        cw = mesh_.lookupObject<volScalarField>("cw");
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find cw in the "
            << "database" << exit(FatalError);
    }
*/
	const tmp<volScalarField>& cw = mesh_.lookupObject<volScalarField>("cw");

    word name(type());

    return store(name, calcMassFlux(cw));
}


bool Foam::functionObjects::absorbedMassFluxDensityObject::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& absorbedMassFluxDensity =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& mfp = absorbedMassFluxDensity.boundaryField()[patchi];

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
