/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::InterfaceRegressionPhaseSystem

Description
    Class which models non-thermally-coupled mass transfers; i.e.,
    representation changes, rather than phase changes.

SourceFiles
    InterfaceRegressionPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef InterfaceRegressionPhaseSystem_H
#define InterfaceRegressionPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interfaceTrackingModel;

/*---------------------------------------------------------------------------*\
                    Class InterfaceRegressionPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class InterfaceRegressionPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

    // Protected data

        dictionary phaseProperties;

        // Adiabatic Flame Temperature
        volScalarField Tad;
        volScalarField Hs1;
        volScalarField Hs2;

        // Burning rate
        volScalarField rb_;
        volScalarField mdotProp_;

        // propellant density
        dimensionedScalar rhoPropellant;

        // Regression Control Parameters
        bool regress_;
        word checkAndRegress_;
        scalar stopRegress_;
        scalar intialV_;
        dimensionedScalar MR_;

        // Regresssion Parameters
        volScalarField alphaOld;
        volScalarField regressionAlpha;
        scalar epsilon;

        // Inlet velocity
        bool thinFlimModel_;
        volVectorField Ug_;
        volVectorField Up_;

    // Protected member functions

        //- Return the representation mass transfer rate
        


public:

    // Constructors

        //- Construct from fvMesh
        InterfaceRegressionPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~InterfaceRegressionPhaseSystem();

    typedef PtrListDictionary<phaseModel> phaseModelList;

    // Public Data
    EntrainedSystem eta;

    // Member Functions

    //- Return the saturationModel
    const saturationModel& saturation() const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Regress Propellant Surface
        virtual void solve();

        //- Return propellant mass flow rate
        virtual tmp<volScalarField> mdotProp() const;

        //- Store old Times
        virtual void store();

        //- Correct the mass transfer rates
        virtual void correct();

        //- Calculate velocity
        void calculateVelocity();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InterfaceRegressionPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
