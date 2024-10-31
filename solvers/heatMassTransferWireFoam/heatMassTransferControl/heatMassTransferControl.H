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

Class
    Foam::heatMassTransferControl

Description
	Solution control class for heat mass transfer solving in 
	heatMassTransferPimpleFoam solver.

SourceFiles
    heatMassTransferControlI.H
    heatMassTransferControl.C

\*---------------------------------------------------------------------------*/

#ifndef heatMassTransferControl_H
#define heatMassTransferControl_H

#include "multiRegionSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatMassTransferControl Declaration
\*---------------------------------------------------------------------------*/

class heatMassTransferControl
:
    public multiRegionSolutionControl
{
protected:

    // Protected data

        //- Maximum number of heat mass transfer iterations
        label nHeatMassIter_;

        //- Current heat mass transfer iteration
        label heatMassIter_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("heatMassTransferControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        heatMassTransferControl
		(
			const Time& time, 
			const word& algorithmName = "heatMassTransfer"
		);


    //- Destructor
    virtual ~heatMassTransferControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();


        // Access

            //- Maximum number of heat and mass transfer iterations
            inline label nHeatMassIter() const;

            //- Current number of heat and mass transfer iterations
            inline label heatMassIter() const;

            //- Flag to indicate the first heat and mass transfer iteration
            inline bool firstHeatMassIter() const;

            //- Flag to indicate the last heat and mass transfer iteration
            inline bool finalHeatMassIter() const;


		// Modify

            //- Raise current heat and mass transfer iteration number
            inline bool nextIter();

            //- Reset current heat and mass transfer iteration number
            inline bool resetIter();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "heatMassTransferControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
