label purePropellantSize = 0;
{
  const volScalarField& propellant = phases[propellantIndex];

  forAll(propellant, i)
  {
    if (propellant[i] == 1)
    {
      purePropellantSize++;
    }
  }
}
labelList purePropellantCells(purePropellantSize);
scalar Tad = fluid.get<scalar>("Tad");
scalarField setTemp(purePropellantSize, Tad);
{
  const volScalarField& propellant = phases[propellantIndex];

  label j = 0;
  forAll(propellant, i)
  {
    if (propellant[i] == 1)
    {
      purePropellantCells[j] = i;
      j++;
    }
  }
}
