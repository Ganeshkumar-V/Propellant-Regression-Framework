for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    EEqnProfile.heatTransfer.start();
    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();
    EEqnProfile.heatTransfer.stop();

    forAll(fluid.anisothermalPhases(), anisothermalPhasei)
    {
      EEqnProfile.Construct.start();
        phaseModel& phase = fluid.anisothermalPhases()[anisothermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        const volVectorField& U = phase.U();

        // Check for presence of phase
        if (max(alpha).value() < 1e-11)
        {
          Warning << phase.name() << " is not present"
                  << " Max alpha: " << max(alpha).value() << endl;
          continue;
        }

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvOptions(alpha, rho, phase.thermoRef().he())
        );
      EEqnProfile.Construct.stop();

      EEqnProfile.Relax.start();
        EEqn.relax();
        fvOptions.constrain(EEqn);
      EEqnProfile.Relax.stop();


        // set Propellant phase Values
        // EEqn.setValues
        // (
        //   purePropellantCells,
        //   phase.thermo().he(setTemp, setTemp, purePropellantCells)
        // );
      EEqnProfile.Solve.start();
        EEqn.solve();
        fvOptions.correct(phase.thermoRef().he());
      EEqnProfile.Solve.stop();
    }

  EEqnProfile.correctThermo.start();
    fluid.correctThermo();
  EEqnProfile.correctThermo.stop();

  EEqnProfile.correct.start();
    fluid.correct();
  EEqnProfile.correct.stop();
}
