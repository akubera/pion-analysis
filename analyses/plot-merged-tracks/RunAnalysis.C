
#include <string>

void
RunAnalysis()
{
  auto *input_files = new TChain("aodTree");
  input_files->Add("/alice/data/2015/LHC15o/000246808/pass1/AOD194/0001/AliAOD.root");

  auto mgr = new AliAnalysisManager("mgr");

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/AddTaskFemto.C");
  // macro_config_path = (Path(__file__).parent / "ConfigFemtoAnalysis.C").absolute()

  std::string path = "/home/akubera/Physics/pion-analysis/analyses/plot-merged-tracks/ConfigFemtoAnalysis.C";
  std::string line = std::string("AddTaskFemto(\"") + path + std::string("\")");

  std::cout << "line: " << line << "\n";
  int err = 0;
  gROOT->ProcessLine(line.c_str(), &err);
  std::cerr << "Error: " << err << "\n";
  // AddTaskFemto(path);

  mgr->InitAnalysis();
  mgr->Print();
  mgr->StartAnalysis("local", input_files);
}
