using Process = std::vector <Insertion>;
using Processes = std::vector <Process>;

// This function takes as input a specific process, and returns a std::string, 
// that corresponds to a name for this process.
// This function can serve many purposes. In our code we use it to assign name  
// to the functions in the numerical library.
std::string processName(Process const &proc)
{
   std::ostringstream sout;
   for (const auto &ins : proc)
         if (ins.isIncoming()) 
         {
            if (!ins.isParticle()) sout << "anti_";
            sout << ins.getField()->getName();
            sout << '_';
         }
   sout << "to_";
   for (const auto &ins : proc)
         if (!ins.isIncoming()) 
         {
            if (!ins.isParticle()) sout << "anti_";
            sout << ins.getField()->getName();
            sout << '_';
         }
   std::string name = sout.str();
   name.erase(name.end() - 1); // erasing the last '-'
   return name;
}

typedef struct
{
  Process process;
  mty::Order order;         // Options mty::Order::TreeLevel, mty::Order::OneLoop, 
  bool leading_order;       // true N1 N1 -> gamma gamma
  mty::gauge::Type Wgauge;   
  
  void printname()
  {
      std::cout << processName(process) << std::endl;
  }

}Process2to2ToCompute;
