///
/// \file mrc-3dratio/RunAnalysis.C
///

#include "../analysis-runner.hh"

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


void
setup_grid(AliAnalysisManager *mgr, TString workdir, TString grid_mode)
{
  auto *alien = new AliAnalysisAlien();
  alien->SetRunMode(grid_mode);
  alien->SetGridOutputDir("output");
  alien->SetGridWorkingDir(workdir);
  alien->SetAliPhysicsVersion("vAN-20190808_ROOT6-1");
  alien->SetDropToShell(false);
  alien->SetCheckCopy(false);
  alien->SetMaxMergeFiles(20);
  alien->SetMaxMergeStages(10);
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
  alien->AddDataFile("/alice/cern.ch/user/a/akubera/xml/4f764ae3df84662c7af07f75f3cb4a71.xml");

  mgr->SetGridHandler(alien);
}


void
usage(std::ostream &out)
{
  out << "Usage:\n";
}


void
AddTasks()
{
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE, kTRUE, kTRUE)");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/AddNuTaskPionPionRoot6.C+");

  gROOT->ProcessLine(R""(
    AddNuTaskPionPionRoot6(
      "true_pions",
      "macro='%%/ConfigNuFemtoAnalysisR6.C+'",
      R"(
         +m; {0:90}; (0.2:0.3,0.6:0.7,0.8:1.0);
         ~do_kt_trueq3d_cf = false; ~do_trueqinv_cf = false; ~do_kt_trueqinv_cf = true;
         ~do_sharequality_cf=false; ~q3d_bin_count=39; ~q3d_maxq=0.117;
         ~mcwg_lednicky = true; ~mcwg_strong = false; ~mcwg_3body = false;
         ~eventreader_filter_bit=128;  ~eventreader_multibit=true;  ~eventreader_read_full_mc=false;
       ~eventreader_dca_globaltrack = 0;
        $pion_1_max_impact_xy = 2.0;
        $pion_1_max_impact_z = 2.0;
         @is_mc_analysis = true; ~trueq3d_extra_bins = true; @num_events_to_mix = 1;
         @enable_pair_monitors=false;  $mc_pion_only=true;
         $pion_1_ideal_pid=211;
         $pair_delta_eta_min = 0.01; $pair_delta_phi_min = 0.03;
         $pion_1_rm_neg_lbl = false; $pion_1_sigma = 3.0;
         $pion_1_min_tpc_chi_ndof = 0.33; )"); )"");

  gROOT->ProcessLine(R""(
    AddNuTaskPionPionRoot6(
      "true_kaons",
      "macro='%%/ConfigNuFemtoAnalysisR6.C+'",
      R"(
         +m; {0:90}; (0.2:0.3,0.6:0.7,0.8:1.0);
         ~do_kt_trueq3d_cf = false; ~do_trueqinv_cf = false; ~do_kt_trueqinv_cf = true;
         ~do_sharequality_cf=false; ~q3d_bin_count=39; ~q3d_maxq=0.117;
         ~mcwg_lednicky = true; ~mcwg_strong = false; ~mcwg_3body = false;
         ~eventreader_filter_bit=128;  ~eventreader_multibit=true;  ~eventreader_read_full_mc=false;
       ~eventreader_dca_globaltrack = 0;
        $pion_1_max_impact_xy = 2.0;
        $pion_1_max_impact_z = 2.0;
         @is_mc_analysis = true; ~trueq3d_extra_bins = true; @num_events_to_mix = 1;
         @enable_pair_monitors=false;  $mc_pion_only=true;
         $pion_1_ideal_pid=321;
         $pair_delta_eta_min = 0.01; $pair_delta_phi_min = 0.03;
         $pion_1_rm_neg_lbl = false; $pion_1_sigma = 3.0;
         $pion_1_min_tpc_chi_ndof = 0.33; )"); )"");


  return;

  gROOT->ProcessLine(R""(
  AddNuTaskPionPionRoot6(
    "mrcvary_fb128_2",
    "macro='%%/ConfigNuFemtoAnalysisR6.C+'",
    R"(
       +m; {0:90}; (0.2:0.3:0.4,0.6:0.7,0.8:1.0);
       ~do_kt_trueq3d_cf = true; ~do_trueqinv_cf = true; ~do_kt_trueqinv_cf = true;
       ~do_sharequality_cf=true; ~q3d_bin_count=39; ~q3d_maxq=0.117;
       ~mcwg_lednicky = true; ~mcwg_strong = false; ~mcwg_3body = false;
       ~eventreader_filter_bit=128;  ~eventreader_multibit=true;  ~eventreader_read_full_mc=true;
       @is_mc_analysis = true; ~trueq3d_extra_bins = true; @num_events_to_mix = 4;
       @enable_pair_monitors=false;  $mc_pion_only=true; $pair_delta_eta_min = 0.01;
       $pair_delta_phi_min = 0.03; $pion_1_sigma = 3.0; $pion_1_min_tpc_chi_ndof = 0.33;
       $pion_1_rm_neg_lbl = true; )"); )"");

  gROOT->ProcessLine(R""(
  AddNuTaskPionPionRoot6(
    "mrcvary_fb128_3",
    "macro='%%/ConfigNuFemtoAnalysisR6.C+'",
    R"(
       +m; {0:90}; (0.2:0.3:0.4,0.6:0.7,0.8:1.0);
       ~do_kt_trueq3d_cf = true; ~do_trueqinv_cf = true; ~do_kt_trueqinv_cf = true;
       ~do_sharequality_cf=true; ~q3d_bin_count=39; ~q3d_maxq=0.117;
       ~mcwg_lednicky = true; ~mcwg_strong = false; ~mcwg_3body = false;
       ~eventreader_filter_bit=128;  ~eventreader_multibit=true;  ~eventreader_read_full_mc=true;
       ~eventreader_dca_globaltrack = 0;
       @is_mc_analysis = true; ~trueq3d_extra_bins = true; @num_events_to_mix = 4;
       @enable_pair_monitors=false;  $mc_pion_only=true; $pair_delta_eta_min = 0.01;
       $pair_delta_phi_min = 0.03; $pion_1_sigma = 3.0; $pion_1_min_tpc_chi_ndof = 0.33;
       $pion_1_max_impact_xy=0.02; $pion_1_max_impact_z=0.04;
       $pion_1_rm_neg_lbl = true; )"); )"");

}


void
RunLocal(TString wd)
{
  TStopwatch timer;
  timer.Start();

  auto *mgr = NewLocalAnalysisManager();

  const TString
    timestamp = wd(wd.Index('-')+1, 14),
    output_filename = Form("MrcResult-%s.root", timestamp.Data());

  mgr->SetCommonFileName(output_filename);

  gSystem->mkdir(wd);
  gSystem->CopyFile("ConfigFemtoAnalysis.C", wd + "/ConfigFemtoAnalysis.C");
  gSystem->CopyFile("RunAnalysis.C", wd + "/RunAnalysis.C");
  gSystem->cd(wd);

  AddTasks();

  mgr->InitAnalysis();
  mgr->PrintStatus();

  TChain *input = new TChain("aodTree");
  // for (int run_num : {33, 34, 35, 36, 67, }) {
    // input->Add(Form("/alice/sim/2016/LHC16g1a/246980/AOD/%03d/AliAOD.root", run_num));
    // input->Add(Form("/alice/sim/2016/LHC16i3a/246982/AOD/%03d/AliAOD.root", run_num));
  for (int run_num : {1, 2, 4, 5, 7,8,9}) {
    input->Add(Form("/alice/sim/2016/LHC16i3a/246982/AOD198/%04d/AliAOD.root", run_num));
  }

  mgr->StartAnalysis("local", input);
  timer.Stop();
  timer.Print();

  TString outfile = wd + "/" + mgr->GetCommonFileName();
  std::cout << "Output written to " << outfile << "\n";
}


void
RunGrid(TString wd)
{
  auto *grid = TGrid::Connect("alien://");
  if (!grid) {
    std::cerr << "Could not connect to alien\n";
    return;
  }

  auto *mgr = NewGridRunAnalysisManager(wd);

  setup_grid(mgr, wd, "full");

  mgr->InitAnalysis();
  mgr->PrintStatus();

  mgr->StartAnalysis("grid");
}


void
RunMerge(TString wd)
{
  auto *grid = TGrid::Connect("alien://");
  if (!grid) {
    std::cerr << "Could not connect to alien\n";
    return;
  }

  auto *mgr = NewGridMergeAnalysisManager(wd);

  setup_grid(mgr, wd, "terminate");

  mgr->InitAnalysis();

  mgr->StartAnalysis("merge");
}


void
// RunAnalysis(TString wd="")
RunAnalysis(std::vector<std::string> args)
{
  if (args.size() == 0) {
    usage(std::cerr);
    return;
  }

  if (args[0] == "local") {
    TString wd = args.size() > 1
               ? TString(args[1])
               : get_timestamp("local-");

    RunLocal(wd);
    return;
  }

  if (args[0] == "grid") {

  }

  if (args[0] == "merge") {
    if (args.size() == 1) {
      usage(std::cerr);
      return;
    }

    TString wd = args[1];
    RunMerge(wd);
    return;
  }

  // TString output_filename = Form("MrcResult-%s.root", timestamp.Data());
}
