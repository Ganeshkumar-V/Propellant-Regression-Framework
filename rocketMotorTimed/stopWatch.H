#include <ctime>

class stopWatch
{
	private:

      std::clock_t startTime;
		  double totalTime;

	public:

      stopWatch()
		  {
		  	   totalTime = 0;
		  }
		  void start()
		  {
		  	   startTime = std::clock();
		  }
		  void stop()
		  {
		  	   totalTime += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
		  }
		  double getTotalTime()
		  {
		  	   return totalTime;
		  }
};
