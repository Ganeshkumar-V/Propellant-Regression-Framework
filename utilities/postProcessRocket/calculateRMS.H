const scalarField UgSqr(UgasF&UgasF);
const scalarField rhogSqr(sqr(rhoGas));
const scalarField UpSqr(UparticlesF&UparticlesF);
const scalarField rhopSqr(sqr(rhoParticles));

const scalarField UgasFz(UgasF.component(vector::Z));
const scalarField UparticlesFz(UparticlesF.component(vector::Z));

const tmp<volScalarField> tgammag(phases[0].thermo().gamma());
const volScalarField& Gammag(tgammag());
const scalarField gammagF(Gammag.boundaryField()[outIndex]);

volScalarField nutGas
(
    IOobject
    (
      "nut.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);
volScalarField nutParticles
(
    IOobject
    (
      "nut.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volScalarField kGas
(
    IOobject
    (
      "k.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);
volScalarField ThetaParticles
(
    IOobject
    (
      "Theta.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

const scalarField& nutgF(nutGas.boundaryField()[outIndex]);
const scalarField& nutpF(nutParticles.boundaryField()[outIndex]);

const scalarField& kgF(kGas.boundaryField()[outIndex]);
const scalarField& ThetapF(ThetaParticles.boundaryField()[outIndex]);

const volScalarField divUg(fvc::div(Ugas));
const scalarField divUgF(divUg.boundaryField()[outIndex]);
const volScalarField divUp(fvc::div(Uparticles));
const scalarField divUpF(divUp.boundaryField()[outIndex]);

tmp<volTensorField> tgradUg = fvc::grad(Ugas);
const volScalarField gradUgzz(tgradUg().component(tensor::ZZ));
const scalarField gradUgzzF(gradUgzz.boundaryField()[outIndex]);

tmp<volTensorField> tgradUp = fvc::grad(Uparticles);
const volScalarField gradUpzz(tgradUp().component(tensor::ZZ));
const scalarField gradUpzzF(gradUpzz.boundaryField()[outIndex]);

const scalarField UgPrimeSqr((2.0/3.0)*(kgF/alphaGasF) - 2.0*(nutgF/alphaGasF)*(gradUgzzF - (1.0/3.0)*divUgF));
const scalarField UpPrimeSqr(ThetapF);

const scalarField ReygF(UgPrimeSqr*rhoGas);
const scalarField ReypF(ThetapF*rhoParticles);
Info << "ReygF: " << ReygF[15] << endl;
Info << "UgPrimeSqr: " << UgPrimeSqr[15] << endl;
Info << "kgF: " << kgF[15] << endl;

const scalarField MachgasF(Mach.boundaryField()[outIndex]);
const scalarField rhogPrimeSqr(rhogSqr*sqr((gammagF-1.0)*sqr(MachgasF))*UgPrimeSqr/UgSqr);
const scalarField rhogPrimeUgPrime(rhoGas*(gammagF-1.0)*sqr(MachgasF)*UgPrimeSqr/UgasFz);
const scalarField rhopPrimeSqr(0.0*rhogPrimeSqr);
const scalarField rhopPrimeUpPrime(0.0*rhogPrimeUgPrime);

Info << "rhogPrimeSqr: " << rhogPrimeSqr[15] << endl;
Info << "gammagF: " << gammagF[15] << endl;
Info << "MachgasF: " << MachgasF[15] << endl;
Info << "UgSqr: " << UgSqr[15] << endl;
Info << "UgasFz: " << UgasFz[15] << endl;
Info << "rhogPrimeUgPrime: " << rhogPrimeUgPrime[15] << endl;

// Intermittency of particle phase will be very high if the volume fraction is very low...
// Intermittency should be calculated differently... One assumption is take intermittency same as the gas phase. 
// Since, intermittency is same, cross correlation is effectively zero. 
const scalarField RMSmassFlow
(
    sqrt
    (
        alphaGasF*(1.0 - alphaGasF)*rhogSqr*UgSqr 
        + alphaGasF*rhogPrimeSqr*UgSqr + alphaGasF*rhoGas*ReygF
        + 2.0*alphaGasF*rhoGas*UgasFz*rhogPrimeUgPrime 
        + 1e-9*alphaParticlesF*(1.0 - alphaParticlesF)*rhopSqr*UpSqr 
        + alphaParticlesF*rhopPrimeSqr*UpSqr + alphaParticlesF*rhoParticles*ReypF
        + 2.0*alphaParticlesF*rhoParticles*UparticlesFz*rhopPrimeUpPrime
        - Foam::sqrt(1e-9)*2.0*alphaGasF*rhoGas*UgasFz*alphaParticlesF*rhoParticles*UparticlesFz
    )*mag(Sf)
);
Info << "RMSmassFlow: " << RMSmassFlow[15] << endl;

const volVectorField gradkg(fvc::grad(kGas));
const scalarField gradkgF(gradkg.component(vector::Z)().boundaryField()[outIndex]);
const scalarField UgPrimePPrime(-(2.0/(3.0*0.9))*rhoGas*nutgF*gradkgF);
const scalarField UpPrimePPrime(0.0*UgPrimePPrime);

const scalarField rhogPrimePPrime(rhoGas*(gammagF-1.0)*sqr(MachgasF)*UgPrimePPrime/UgasFz);
const scalarField rhopPrimePPrime(0.0*rhogPrimePPrime);

const scalarField pPrimeSqr(0.34*rhogSqr*sqr(kgF));

Info << "rhogPrimePPrime: " << rhogPrimePPrime[15] << endl;
Info << "UgPrimePPrime: " << UgPrimePPrime[15] << endl;
Info << "pPrimeSqr: " << pPrimeSqr[15] << endl;
Info << "nutgF: " << nutgF[15] << endl;
Info << "gradkgF: " << gradkgF[15] << endl;

const scalarField RMSThrust
(
    sqrt
    (
        alphaGasF*(1.0 - alphaGasF)*rhogSqr*UgSqr*UgSqr + alphaGasF*rhogPrimeSqr*UgSqr*UgSqr
        - alphaGasF*sqr(ReygF) - alphaParticlesF*sqr(ReypF) 
        - Foam::sqrt(1e-9)*2.0*alphaGasF*rhoGas*UgSqr*alphaParticlesF*rhoParticles*UpSqr
        + 1e-9*alphaParticlesF*(1.0 - alphaParticlesF)*rhopSqr*UpSqr*UpSqr 
        + alphaParticlesF*rhopPrimeSqr*UpSqr*UpSqr
        - 2.0*alphaGasF*rhoGas*UgSqr*alphaParticlesF*ReypF - Foam::sqrt(1e-9)*2.0*alphaGasF*ReygF*alphaParticlesF*rhoParticles*UpSqr
        - 2.0*alphaGasF*ReygF*alphaParticlesF*ReypF 
        + pPrimeSqr 
        + 4.0*alphaGasF*rhoGas*UgasFz*UgPrimePPrime + 2.0*alphaGasF*UgSqr*rhogPrimePPrime
        + 4.0*alphaParticlesF*rhoParticles*UparticlesFz*UpPrimePPrime + 2.0*alphaParticlesF*UpSqr*rhopPrimePPrime
        + alphaGasF*rhoGas*UgSqr*(UgasFz*rhogPrimeUgPrime - 2.0*alphaGasF*ReygF + 6.0*ReygF)
        + alphaParticlesF*rhoParticles*UpSqr*(UparticlesFz*rhopPrimeUpPrime - 2.0*alphaParticlesF*ReypF + 6.0*ReypF)
    )*mag(Sf)
);
Info << "RMSThrust: " << RMSThrust[15] << endl;

tFgas += sum(alphaGasF*ReygF*mag(Sf));
tFparticles += sum(alphaParticlesF*ReypF*mag(Sf)); 

const scalarField FprimeMdotprime
(
    (
        alphaGasF*(1.0 - alphaGasF)*rhoGas*UgSqr*UgasFz + alphaGasF*rhogPrimeSqr*UgSqr*UgasFz 
        + alphaGasF*(3.0 - alphaGasF)*rhoGas*UgasFz*ReygF
        + 1e-9*alphaParticlesF*(1.0 - alphaParticlesF)*rhoParticles*UpSqr*UparticlesFz  
        + alphaParticlesF*rhopPrimeSqr*UpSqr*UparticlesFz 
        + 1e-9*alphaParticlesF*(3.0 - alphaGasF)*rhoParticles*UparticlesFz*ReypF 
        - Foam::sqrt(1e-9)*alphaGasF*rhoGas*UgasFz*alphaParticlesF*rhoParticles*UparticlesFz*(UgasFz + UparticlesFz)
        - alphaGasF*alphaParticlesF*
            (
                rhoGas*UgasFz*ReypF + rhoParticles*UparticlesFz*ReygF 
                + 3.0*alphaGasF*rhoGas*UgSqr*rhogPrimeUgPrime + 3.0*alphaParticlesF*rhoParticles*UpSqr*rhopPrimeUpPrime
            )
        + alphaGasF*rhoGas*UgPrimePPrime + alphaGasF*UgasFz*rhogPrimePPrime 
        + alphaParticlesF*rhoParticles*UpPrimePPrime + alphaParticlesF*UparticlesFz*rhopPrimePPrime
    )*sqr(mag(Sf))
);
Info << "FprimeMdotprime: " << FprimeMdotprime[15] << endl;

tmdotRMS = sum(RMSmassFlow);
tFRMS = sum(RMSThrust);
tFprimeMdotprime = sum(FprimeMdotprime);

// Isp and Isp RMS
tIspRMS = tIsp * Foam::sqrt(sqr(tFRMS/tF) + sqr(tmdotRMS/tmdot) - 2.0*tFprimeMdotprime/(tF*tmdot));