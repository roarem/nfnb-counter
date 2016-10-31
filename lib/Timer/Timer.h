#pragma once
#include <chrono>

class Timer
{
  public:

    Timer (){};

    void    startTimer	      ()  {my_start = clock::now();}
    void    stopTimer	      ()  {my_end = clock::now();}
    //int     elapsedTimeSeconds();
    //int     elapsedTimeMilli  ();
    //int     elapsedTimeMicro  ();

    int elapsedTimeSeconds()
    {
      int count = std::chrono::duration_cast<seconds_type>(my_end - my_start).count();
      return count;
    }
    
    int elapsedTimeMilli()
    {
      int count = std::chrono::duration_cast<milli_type>(my_end - my_start).count();
      return count;
    }
    
    int elapsedTimeMicro()
    {
      int count = std::chrono::duration_cast<micro_type>(my_end - my_start).count();
      return count;
    }

    int *elapsedTimeClock()
    {
        int seconds = elapsedTimeSeconds();
        int *returnTime = new int [3];
        returnTime[2] = seconds%60;
        returnTime[1] = seconds/60%60;
        returnTime[0] = seconds/60/60%60;
        return returnTime; 
    }

  protected:

    class System* my_system = nullptr;

    typedef std::chrono::high_resolution_clock	clock;
    typedef std::chrono::seconds        seconds_type;
    typedef std::chrono::milliseconds   milli_type;
    typedef std::chrono::microseconds   micro_type;

    std::chrono::duration<int> seconds;
    std::chrono::duration<int> milli;
    std::chrono::duration<int> micro;

    clock::time_point	      my_start;
    clock::time_point	      my_end;
};



