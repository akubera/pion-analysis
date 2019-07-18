///
/// \file mrc-hypercube/Analysis.C
///

#include <iostream>
#include <string>

/// \class Analysis
///
struct Analysis {
  // const std::string name {"MrcHypercube"};
  static std::string name()
    { return "MrcHypercube"; }

  // std::()
};


void SetupAnalysis()
{
  std::cout << "[SetupAnalysis]\n";
}
