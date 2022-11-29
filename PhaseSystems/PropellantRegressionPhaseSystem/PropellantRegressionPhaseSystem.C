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

\*---------------------------------------------------------------------------*/

#include "PropellantRegressionPhaseSystem.H"
#include "interfaceTrackingModel.H"
#include "fvmSup.H"
#include "phaseSystem.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::rDmdt
(
    const phasePairKey& key
) const
{
    if (!rDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar rDmdtSign(Pair<word>::compare(rDmdt_.find(key).key(), key));
    return rDmdtSign**rDmdt_[key];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::PropellantRegressionPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    )
{
    this->generatePairsAndSubModels
    (
        "interfaceTracking",
        interfaceTrackingModels_
    );

    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        this->rDmdt_.set
        (
            interfaceTrackingModelIter.key(),
            phaseSystem::dmdt(interfaceTrackingModelIter.key()).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::
~PropellantRegressionPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return BasePhaseSystem::dmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    // Fill the mass transfer rates with zero
    this->fillFields("dmdt", dimDensity/dimTime, dmdts);

    forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();
        // Add massTransfer Rate to the Gas Phase only
        this->addField(pair.phase2(), "dmdt", rDmdt, dmdts);
    }
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.set
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

    // Mass transfer across the interface
    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[interfaceTrackingModelIter.key()]);

        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        // Note that the phase YiEqn does not contain a continuity error term,
        // so these additions represent the entire mass transfer

        const volScalarField dmdt(this->rDmdt(pair));
        const volScalarField dmdt12(negPart(dmdt));
        const volScalarField dmdt21(posPart(dmdt));

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            const word name
            (
                IOobject::groupName(Yi[i].member(), phase.name())
            );

            const word otherName
            (
                IOobject::groupName(Yi[i].member(), otherPhase.name())
            );

            *eqns[name] +=
                dmdt21*eqns[otherName]->psi()
              + fvm::Sp(dmdt12, eqns[name]->psi());

            *eqns[otherName] -=
                dmdt12*eqns[name]->psi()
              + fvm::Sp(dmdt21, eqns[otherName]->psi());
        }

    }

    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
  autoPtr<phaseSystem::heatTransferTable> eqnsPtr
  (
      new phaseSystem::heatTransferTable()
  );

  phaseSystem::heatTransferTable& eqns = eqnsPtr();

  forAll(this->phaseModels_, phasei)
  {
      const phaseModel& phase = this->phaseModels_[phasei];

      eqns.set
      (
          phase.name(),
          new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
      );
  }

  forAllConstIter
  (
      interfaceTrackingModelTable,
      interfaceTrackingModels_,
      interfaceTrackingModelIter
  )
  {
      // const volScalarField qg(interfaceTrackingModelIter()->Qg());1
      // const volScalarField ql(interfaceTrackingModelIter()->Ql());

      const phasePair& pair(this->phasePairs_[interfaceTrackingModelIter.key()]);

      if (pair.ordered())
      {
          continue;
      }

      const phaseModel& phase1 = pair.phase1();
      const phaseModel& phase2 = pair.phase2();

      const dimensionedScalar one(dimless, 1.0);

      // dmdt term is accounted in the continuity error term
      // (*eqns[phase2.name()]).source() += qg;
      // (*eqns[phase1.name()]).source() -= ql;

      *eqns[phase2.name()] += fvm::laplacian
      (
          fvc::interpolate(phase2.alphaEff()),
          phase2.thermo().he()
      ) - fvm::laplacian
      (
        fvc::interpolate(phase2)*fvc::interpolate(phase2.alphaEff()),
        phase2.thermo().he()
      );

      *eqns[phase1.name()] += fvm::laplacian
      (
          fvc::interpolate(phase1.alphaEff()),
          phase1.thermo().he()
      ) - fvm::laplacian
      (
        fvc::interpolate(phase1)*fvc::interpolate(phase1.alphaEff()),
        phase1.thermo().he()
      );

  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::solve()
{
  // Regress Propellant surface (Manipulate propellant volume fraction)
  forAllIter
  (
    interfaceTrackingModelTable,
    interfaceTrackingModels_,
    interfaceTrackingModelIter
  )
  {
    word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
    volScalarField& alpha = this->db().template lookupObjectRef<volScalarField>(propellant);
    interfaceTrackingModelIter()->regress(alpha);
  }

  // Solve other phase volume fraction equations
  BasePhaseSystem::solve();
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        *rDmdt_[interfaceTrackingModelIter.key()] =
            dimensionedScalar(dimDensity/dimTime);
    }

    //- Finds burning rate (rb = aP^n)
    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        interfaceTrackingModelIter()->correct();
    }

    //- return burning Rate
    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
      interfaceTrackingModelIter()->rb();
    }

    //- find mass source terms

}

template<class BasePhaseSystem>
bool Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
