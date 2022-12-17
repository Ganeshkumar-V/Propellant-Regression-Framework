#include "stopWatch.H"

class EEqntimeProfiler
{
  public:
    stopWatch heatTransfer;
    stopWatch Construct;
    stopWatch Relax;
    stopWatch Solve;
    stopWatch correctThermo;
    stopWatch correct;
    void displayTime()
    {
      Info << "EEqn Time Profile: " << endl;
      Info << "     ->  heatTransfer : " << heatTransfer.getTotalTime() << endl;
      Info << "     ->     Construct : " << Construct.getTotalTime() << endl;
      Info << "     ->         Relax : " << Relax.getTotalTime() << endl;
      Info << "     ->         Solve : " << Solve.getTotalTime() << endl;
      Info << "     -> correctThermo : " << correctThermo.getTotalTime() << endl;
      Info << "     ->       correct : " << correct.getTotalTime() << endl;
      Info << "     ->         Total = " << correct.getTotalTime()
                                            + correctThermo.getTotalTime()
                                            + Solve.getTotalTime()
                                            + Relax.getTotalTime()
                                            + Construct.getTotalTime()
                                            + heatTransfer.getTotalTime() << endl;
    }
};
