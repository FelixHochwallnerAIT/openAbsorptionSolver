/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "pimpleMultiFluidRegionControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleMultiFluidRegionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::pimpleMultiFluidRegionControl::time
(
    const PtrList<fvMesh>& pimpleMeshes
)
{
    if (pimpleMeshes.empty())
    {
        FatalErrorInFunction
            << "There needs to be at least one region"
            << exit(FatalError);
    }

    return pimpleMeshes[0].time();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleMultiFluidRegionControl::pimpleMultiFluidRegionControl
(
    PtrList<fvMesh>& pimpleMeshes,
    PtrList<word>& solveRegionFluidFlow,
    const word& algorithmName
)
:
    multiRegionSolutionControl(time(pimpleMeshes), algorithmName),
    pimpleLoop(static_cast<solutionControl&>(*this)),
    convergenceControl(static_cast<solutionControl&>(*this)),
    correctorConvergenceControl
    (
        static_cast<solutionControl&>(*this),
        "outerCorrector"
    ),
    pimpleControls_(),
    solveRegion(solveRegionFluidFlow)
{
    bool allSteady = true, allTransient = true;

    forAll(pimpleMeshes, i)
    {
        pimpleControls_.append
        (
            new pimpleNoLoopControl(pimpleMeshes[i], algorithmName, *this)
        );

        allSteady = allSteady && pimpleMeshes[i].steady();
        allTransient = allTransient && pimpleMeshes[i].transient();
    }

    read();

    forAll(pimpleMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
        pimpleControls_[i].printResidualControls();

        if (nCorrPimple_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
            pimpleControls_[i].printCorrResidualControls(nCorrPimple_);
        }
    }

    Info<< nl << algorithmName << ": Operating solver in "
        << (allSteady ? "steady-state" : allTransient ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorrPimple_
        << " outer corrector" << (nCorrPimple_ == 1 ? "" : "s") << nl;

    if ((allSteady || allTransient) && nCorrPimple_ == 1)
    {
        Info<< algorithmName << ": Operating solver in "
            << (allSteady ? "SIMPLE" : "PISO") << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleMultiFluidRegionControl::~pimpleMultiFluidRegionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleMultiFluidRegionControl::read()
{
    forAll(pimpleControls_, i)
    {
        if (!pimpleControls_[i].read())
        {
            return false;
        }
    }

    const dictionary& solutionDict = dict();

    nCorrPimple_ = solutionDict.lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


bool Foam::pimpleMultiFluidRegionControl::hasResidualControls() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = result && pimpleControls_[i].hasResidualControls();
    }

    return result;
}


bool Foam::pimpleMultiFluidRegionControl::hasCorrResidualControls() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = result && pimpleControls_[i].hasCorrResidualControls();
    }

    return result;
}


bool Foam::pimpleMultiFluidRegionControl::criteriaSatisfied() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        if(solveRegion[i] == "true")
        {
            result = pimpleControls_[i].criteriaSatisfied() && result;
        }
    }

    return result;
}


bool Foam::pimpleMultiFluidRegionControl::corrCriteriaSatisfied() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        if(solveRegion[i] == "true")
        {
            result = pimpleControls_[i].corrCriteriaSatisfied() && result;
        }
    }

    return result;
}


void Foam::pimpleMultiFluidRegionControl::resetCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].resetCorrSolveIndex();
    }
}


void Foam::pimpleMultiFluidRegionControl::updateCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateCorrSolveIndex();
    }
}


bool Foam::pimpleMultiFluidRegionControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].updateFinal();
        }

        return false;
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].storePrevIterFields();
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateFinal();
    }

    return true;
}


bool Foam::pimpleMultiFluidRegionControl::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].storePrevIterFields();
        }
    }

    return time.run();
}


bool Foam::pimpleMultiFluidRegionControl::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].storePrevIterFields();
        }
    }

    return time.loop();
}


// ************************************************************************* //
