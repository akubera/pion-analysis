///
/// analyses/analysis-runner.hh
///

#pragma once

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TApplication.h>
#include <TGrid.h>
#include <TString.h>
#include <TDatime.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TPython.h>

#include <AliAnalysisAlien.h>
#include <AliAnalysisManager.h>

#include <vector>
#include <string>
#include <iostream>
#include <optional>


// namespace std {
// template <typename T> using optional = std::experimental::optional<T>;
// }

/// Define this in the RunAnalysis macro
void RunAnalysis(std::vector<std::string> args);


void
RunAnalysis()
{
  char **argv = gApplication->Argv();
  const int argc = gApplication->Argc();

  gSystem->SetBuildDir("~/.cache/rootbuild");

  std::vector<std::string> all_args(argv, argv+argc);

  std::string macro_name(gInterpreter->GetCurrentMacroName());
  std::string macro_path(macro_name.begin() + macro_name.find("/./") + 3, macro_name.end());

  std::vector<std::string> args(std::find(all_args.begin(), all_args.end(), macro_path) + 1, all_args.end());
  RunAnalysis(args);
}


AliAnalysisManager*
NewAnalysisManagerExec(TString wd="")
{
  TString abc = wd.IsWhitespace() ? "abc" : wd;
  return new AliAnalysisManager(abc);
}


AliAnalysisManager*
NewAnalysisManagerMerge(TString wd="")
{
  TString abc = wd.IsWhitespace() ? "abc" : wd;
  return new AliAnalysisManager("mergemgr");
}


TString get_timestamp(TString prefix="", TString suffix="")
{
  TDatime td;
  TString timestamp = Form("%08d%06d", td.GetDate(), td.GetTime());
  return prefix + timestamp + suffix;
}
