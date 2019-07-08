///
/// \file mrc-1d/RunAnalysis.C
///

#include <TString.h>
#include <TDatime.h>
#include <TChain.h>
#include <TStopwatch.h>

#include <AliAnalysisAlien.h>
#include <AliAnalysisManager.h>

#include <set>
#include <string>
#include <fstream>
#include <iostream>


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

#define MERGE_MODE true


void
setup_grid(AliAnalysisManager *mgr, TString workdir)
{
  auto *alien = new AliAnalysisAlien();
#if MERGE_MODE
  alien->SetRunMode("terminate");
#else
  alien->SetRunMode("full");
#endif
  alien->SetGridOutputDir("output");
  alien->SetGridWorkingDir(workdir);
  alien->SetAliPhysicsVersion("vAN-20190705_ROOT6-1");
  alien->SetDropToShell(false);
  alien->SetCheckCopy(false);
  alien->SetMaxMergeFiles(17);
  alien->SetMaxMergeStages(3);
  alien->SetSplitMaxInputFileNumber(15);
  alien->SetNrunsPerMaster(30);
  alien->SetMergeViaJDL(true);
  alien->SetTTL(12 * 60 * 60);
  alien->AddAdditionalLibrary("ConfigFemtoAnalysis.C");

  // alien->SetGridDataDir("/alice/sim/2016/LHC16i3a");
  // alien->SetGridDataDir("/alice/sim/2018/LHC18e1");
  // alien->SetDataPattern("/AOD/*/AliAOD.root");

  // std::set<int> runs = {
  // 246675, 246676, 246805, 246804, 246807, 246424, 246808,
  // 246809, 246428, 246431, 246945, 246434,
  // 246948, 246844, 246845, 246846, 246847, 246851, 246980, 246982, 246984, 246989, 246991, 246865,
  // 246994, 246867, 246487, 246871, 246488, 246493, 246750, 246751, 246495, 246757, 246758, 246759,
  // 246760, 246763, 246765
  // };

  // std::set<int> negfield = {
  //   246675, 246676, 246805, 246804, 246807, 246424, 246808, 246809, 246428, 246431, 246945, 246434,
  //   246948, 246844, 246845, 246846, 246847, 246851, 246980, 246982, 246984, 246989, 246991, 246865,
  //   246994, 246867, 246487, 246871, 246488, 246493, 246750, 246751, 246495, 246757, 246758, 246759,
  //   246760, 246763, 246765
  // };

  // std::set<int> posfield = {
  //   246272, 246275, 246148, 246276, 245766, 246151, 246152, 246153, 245775, 246036, 246037, 245785,
  //   246042, 246048, 245793, 246178, 246049, 246053, 246182, 246181, 245683, 245692, 245949, 245952,
  //   245954, 245700, 245829, 246087, 245831, 245705, 246089, 246217, 245833, 245963, 246222, 246225,
  //   246113, 245729, 245731, 246115, 245738, 246001, 246003, 245752, 246012, 245759
  // };

  // for (auto run : runs) {
  //   alien->AddRunNumber(run);
  // }

  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/xml/ae190dbdb6da2271f8dd14c1e5b8536d.xml");
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/job-20190623181107/246675_246808.xml");
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/job-20190630183405/246844_246991.xml");
  alien->AddDataFile("/alice/cern.ch/user/a/akubera/xml/a404dfae1b1e3b9b96e69e59a5fd0d8b.xml");

  mgr->SetGridHandler(alien);
}


void
RunAnalysis(TString wd="")
{
  std::cout << "RunAnalysis\n";

  TStopwatch timer;
  timer.Start();

  TDatime td;
  TString timestamp = Form("%08d%06d", td.GetDate(), td.GetTime());

  auto *mgr = new AliAnalysisManager();

  if (wd.IsWhitespace()) {
    wd = "job-" + timestamp;
  } else {
    timestamp = wd(wd.Index('-')+1, 14);
  }

  // std::cout << timestamp << "\n" << wd << "\n";
  gSystem->mkdir(wd);
  gSystem->CopyFile("ConfigFemtoAnalysis.C", wd + "/ConfigFemtoAnalysis.C");
  gSystem->cd(wd);

  TString output_filename = Form("MrcResult-%s.root", timestamp.Data());
  mgr->SetCommonFileName(output_filename);

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE, kTRUE, kTRUE)");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

  auto *femtotask = new AliAnalysisTaskFemtoNu("Mrc1D", "ConfigFemtoAnalysis.C", "");
  mgr->AddTask(femtotask);
  femtotask->SetupContainers();

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
    input->Add(Form("/alice/sim/2016/LHC16g1a/246980/AOD/%03d/AliAOD.root", run_num));
  }

  mgr->StartAnalysis("local", input);
#endif

  timer.Stop();
  timer.Print();

#if !RUN_GRID
  TString outfile = wd + "/" + mgr->GetCommonFileName();
  std::cout << "Output written to " << outfile << "\n";
#endif
}
