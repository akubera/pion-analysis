///
/// \file RunAnalysis.C
///

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TApplication.h>
#include <TString.h>
#include <TDatime.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TPython.h>

#include <AliAnalysisAlien.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskFemtoNu.h>

#include <set>
#include <string>
#include <fstream>
#include <iostream>

#define ALIPHYSICS_VERSION "vAN-20190803_ROOT6-1"


void usage(std::ostream &out)
{
   out <<    "usage: RunAnalysis.C <mode> <options>\n"
       << "\n"
       << "\n    mode:"
       << "\n";
}


std::set<TString>
load_files(std::string filename)
{
  std::set<TString> result;

  std::ifstream f(filename);

  std::string s;
  while (std::getline(f, s)) {
    result.emplace(s.c_str());
  }

  return result;
}


void
setup_grid(AliAnalysisManager *mgr, TString workdir, TString gridmode)
{
  auto *alien = new AliAnalysisAlien();
  alien->SetRunMode(gridmode);
  alien->SetGridOutputDir("output");
  alien->SetGridWorkingDir(workdir);
  alien->SetAliPhysicsVersion(ALIPHYSICS_VERSION);
  alien->SetDropToShell(false);
  alien->SetCheckCopy(false);
  alien->SetMaxMergeFiles(7);
  alien->SetMaxMergeStages(9);
  alien->SetSplitMaxInputFileNumber(15);
  alien->SetNrunsPerMaster(30);
  alien->SetMergeViaJDL(true);
  alien->SetTTL(12 * 60 * 60);
  alien->AddAdditionalLibrary("ConfigFemtoAnalysis.C");

  alien->SetGridDataDir("/alice/sim/2016/LHC16g1");
  alien->SetDataPattern("/AOD198/*/AliAOD.root");

  std::set<int> runs = {
    // negfield
    246945, 246434, 246948, 246844, 246845, 246846, 246847, 246851, 246980, 246982,
    246984, 246989, 246991, 246865, 246994, 246867, 246487, 246871, 246488, 246493,

    // posfield
    // 245683, 245692, 245949, 245952, 245954, 245700, 245829, 246087, 245831, 245705,
    // 246089, 246217, 245833, 245963, 246222, 246225, 246113, 245729, 245731, 246115,
  };

  for (auto run : runs) {
    alien->AddRunNumber(run);
  }
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/xml/ae190dbdb6da2271f8dd14c1e5b8536d.xml");
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/job-20190623181107/246675_246808.xml");
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/job-20190630183405/246844_246991.xml");

  mgr->SetGridHandler(alien);
}



AliAnalysisManager*
CreateRunManager(TString filename)
{
  auto *mgr = new AliAnalysisManager();

  TDatime td;
  TString timestamp = Form("%08d%06d", td.GetDate(), td.GetTime());
  mgr->SetCommonFileName(filename);

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE, kTRUE, kTRUE)");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

  auto *femtotask = new AliAnalysisTaskFemtoNu("femtotask", "ConfigFemtoAnalysis.C", "");
  mgr->AddTask(femtotask);
  femtotask->SetupContainers();

  return mgr;
}


AliAnalysisManager*
CreateMergeManager(TString filename)
{
  auto *mgr = new AliAnalysisManager();
  mgr->SetCommonFileName(filename);

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE, kTRUE, kTRUE)");

  return mgr;
}


TString result_filename_from_wd(const TString wd)
{
  TString timestamp = wd(wd.Index('-')+1, 14);
  TString result = Form("MrcResult-%s.root", timestamp.Data());
  return result;
}

void
RunGridAnalysis(TString wd)
{
  const TString result_filename = result_filename_from_wd(wd);

  auto *mgr = CreateRunManager(result_filename);

  gSystem->mkdir(wd);
  gSystem->CopyFile("ConfigFemtoAnalysis.C", wd + "/ConfigFemtoAnalysis.C");
  gSystem->cd(wd);

  setup_grid(mgr, wd, "full");

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
}


void
RunMergeAnalysis(TString wd)
{
  const TString result_filename = result_filename_from_wd(wd);
  auto *mgr = CreateMergeManager(result_filename);

  setup_grid(mgr, wd, "terminate");

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");

/*
#define RUN_GRID true


#if RUN_GRID
  setup_grid(mgr, wd);

  mgr->InitAnalysis();
  mgr->PrintStatus();

  mgr->StartAnalysis("grid");
#else
  mgr->InitAnalysis();
  mgr->PrintStatus();

  TChain *input = new TChain("aodTree");
  for (int run_num : {1, 2, 3, 4}) {
    input->Add(Form("/alice/sim/2016/LHC16g1/246928/AOD198/%04d/AliAOD.root", run_num));
  }

  mgr->StartAnalysis("local", input);
#endif

  timer.Stop();
  timer.Print();

#if !RUN_GRID
  TString outfile = wd + "/" + mgr->GetCommonFileName();
  std::cout << "Output written to " << outfile << "\n";
#endif
*/
}


void
RunLocalAnalysis(TString wd)
{
}


void
RunAnalysis(std::vector<std::string> args)
{
  if (args.size() < 1) {
    usage(std::cerr);
    return;
  }

  if (args[0] == "merge")  {
    TString wd;

    if (args.size() == 1) {
      TPython::Exec("from pathlib import Path");
      wd = (const char*)TPython::Eval("str(list(Path().glob('job-*'))[-1])");
      std::cout << "Merging latest" << wd <<"\n";
    } else {
      wd = args.at(1);
      std::cout << "Merging: " << wd << "\n";
    }

    RunMergeAnalysis(wd);
    return;
  }

  // create new analysis
  if (args[0] == "run" || args[0] == "new") {

    TDatime td;
    const TString
      timestamp = Form("%08d%06d", td.GetDate(), td.GetTime()),
      workdir = "job-" + timestamp;

    std::cout << "Creating new job: " << workdir << "\n";
    RunGridAnalysis(workdir);
    return;
  }

  // switch (args.size() - 1) {
  // case 1:
  //   RunAnalysis(args[1]);
  //   break;
  // case 2:
  //   // RunAnalysis(args[1], std::stoi(args[2]));
  //   break;
  // default:
    usage(std::cerr);
  // }
}

void
RunAnalysis()
{
  char **argv = gApplication->Argv();
  const int argc = gApplication->Argc();

  std::vector<std::string> all_args(argv, argv+argc);

  std::string macro_name(gInterpreter->GetCurrentMacroName());
  std::string macro_path(macro_name.begin() + macro_name.find("/./") + 3, macro_name.end());

  std::vector<std::string> args(std::find(all_args.begin(), all_args.end(), macro_path) + 1, all_args.end());
  RunAnalysis(args);
}
