/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    heatMassTransferWireFoam

Description
    A multiregion transient heat and mass transfer solver. 
    The fluid flow is of laminar incompressible type.
    The coupled heat and mass transfer problem is solved at each time step,
    after the decoupled flow problem is solved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "icoCourantNo.H"
#include "pimpleMultiFluidRegionControl.H"
#include "heatMassTransferControl.H"
#include "customFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createFluidFields.H"

    #include "initContinuityErrs.H"

    pimpleMultiFluidRegionControl pimples(fluidRegions, solveRegionFluidFlow);

    #include "createTimeControls.H"

    heatMassTransferControl hmtControl(runTime);
    Info<< "Doing " << hmtControl.nHeatMassIter() 
        << " heat and mass transfer iterations." 
        << nl << nl << endl;

    // Read periodic flow field
    IOdictionary periodicFlowDict
    (
        IOobject
        (
            "periodicFlow",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );
    List<word> timesteps = periodicFlowDict.lookup("timesteps");
    PtrList<volVectorField> UPeriodic(timesteps.size());
    Info << "Reading periodic flow field with " << UPeriodic.size() 
        << " timesteps from time " << timesteps[0] 
        << " to " << timesteps[timesteps.size()-1] << "." << nl << endl;    
    forAll(UPeriodic, i)
    {
        Info << "Timestep: " << i+1 
            << ", Reading periodic flow field from time " << timesteps[i] 
            << endl;
        UPeriodic.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    "periodicFlow/" + timesteps[i],
                    fluidRegions[0],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[0]
            )
        );
    }
    Info << "Done reading periodic flow field." << nl << endl;

int i_periodic = 0;
while (pimples.run(runTime))
{
    #include "readTimeControls.H"
    #include "icoMultiRegionCourantNo.H"
    #include "setMultiRegionDeltaT.H"

    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;
    
    // Rereading fluid flow
    forAll(fluidRegions, i)
    {
        if(solveRegionFluidFlow[i] == "true")
        {
            Info<< "\nRereading fluid flow for fluid region "
                << fluidRegions[i].name() << endl;
            #include "rereadPeriodicFlow.H"
        }
        else if(solveRegionFluidFlow[i] == "false")
        {
            Info<< "\nFluid flow in region " 
                << fluidRegions[i].name() 
                << " should not be reread.\n" << endl;
        }
        else
        {
            FatalError  << "Define if fluid flow should be "
                        <<"reread in region (true or false)."
                        << nl << exit(FatalError);
        }
    }

    // Solving heat and mass transfer problem
    Info<< "\n\nSolving heat and mass transfer problem." << endl;
    // Initalize convergence check
    #include "initConvergenceCheckHeatMass.H"
    while (hmtControl.heatMassIter() < hmtControl.nHeatMassIter())
    {
        // Raise iteration number
        hmtControl.nextIter();
        Info<< "\nHeat and mass transfer problem iteration "
                << hmtControl.heatMassIter() << "/" 
                << hmtControl.nHeatMassIter() <<endl;

        forAll(fluidRegions, i)
        {
            Info<< "\nSolving heat and mass transfer problem "
                <<"for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFieldsHeatMass.H"
            #include "readRegionFluidResidualControls.H"
            #include "initRegionFluidConvergenceCheckHeatMass.H"
            #include "solveFluidHeatMass.H"
            #include "saveRegionFluidConvergenceCheckHeatMass.H"
        }
        #include "convergenceCheckHeatMass.H"
    }
    if (hmtControl.nHeatMassIter() == 0)
    {
        Info<< "Heat and mass transfer should not be solved." << endl;
    }
    Info<< nl << endl;
    // Reset iteration number
    hmtControl.resetIter();

    runTime.write();

    // Write fields
    // # include "writeFluidFields.H"
    
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    // Increase periodic flow field
    i_periodic++;
    if (i_periodic == UPeriodic.size())
    {
        i_periodic = 0;
    }
}
Info << "End\n" << endl;

return 0;
}
// ************************************************************************* //
