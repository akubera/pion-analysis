
#include <TString.h>
#include <TDatime.h>
#include <TStopwatch.h>

#include <AliAnalysisManager.h>

#include <iostream>


void
RunAnalysis()
{
  std::cout << "RunAnalysis\n";

  TStopwatch timer;
  timer.Start();

  TDatime td;
  auto *mgr = new AliAnalysisManager();

  TString output_filename = Form("MrcResult-%08d%06d.root", td.GetDate(), td.GetTime());

  mgr->SetCommonFileName(output_filename);

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE, kTRUE, kTRUE)");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");


}
