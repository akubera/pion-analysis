///
/// \file RunAnalysis.C
///

#include <TStopwatch.h>
#include <AliAnalysisManager.h>


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

struct { int x; }["abc"];

void setup_grid(AliAnalysisManager *mgr, TString workdir)
{
  auto *alien = new AliAnalysisAlien();
  alien->SetRunMode("full");
  // alien->SetRunMode("terminate");
  alien->SetGridOutputDir("output");
  alien->SetGridWorkingDir(workdir);
  alien->SetAliPhysicsVersion("vAN-20190605_ROOT6-1");
  alien->SetDropToShell(false);
  alien->SetCheckCopy(false);
  alien->SetMaxMergeFiles(17);
  alien->SetMaxMergeStages(4);
  alien->SetSplitMaxInputFileNumber(23);
  alien->SetNrunsPerMaster(30);
  alien->SetMergeViaJDL(true);
  alien->SetTTL(12 * 60 * 60);
  alien->AddAdditionalLibrary("ConfigFemtoAnalysis.C");

  alien->SetGridDataDir("/alice/sim/2016/LHC16i3a");
  alien->SetDataPattern("/AOD/*/AliAOD.root");

  std::vector<int> runs = {
  246675, 246676, 246805, 246804, 246807, 246424, 246808, 246809, 246428, 246431, 246945, 246434,
  // 246948, 246844, 246845, 246846, 246847, 246851, 246980, 246982, 246984, 246989, 246991, 246865,
  // 246994, 246867, 246487, 246871, 246488, 246493, 246750, 246751, 246495, 246757, 246758, 246759,
  // 246760, 246763, 246765
  };

  for (auto run : runs) {
    alien->AddRunNumber(run);
  }
  // alien->AddDataFile("/alice/cern.ch/user/a/akubera/xml/ae190dbdb6da2271f8dd14c1e5b8536d.xml");

  mgr->SetGridHandler(alien);
}

void
RunAnalysis(TString wd="")
{
  std::cout << "Running Analysis\n";

  TStopwatch timer;
  timer.Start();

  TDatime td;
  auto *mgr = new AliAnalysisManager();

  mgr->SetCommonFileName(Form("QinvResult-%-d%d.root", td.GetDate(), td.GetTime()));
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  // gROOT->Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIUPLICITY/macros/AddTaskMultSelection.C");

  auto *femtotask = new AliAnalysisTaskFemtoNu("femtotask", "ConfigFemtoAnalysis.C", "");
  mgr->AddTask(femtotask);
  femtotask->SetupContainers();

  if (wd.IsWhitespace()) {
    wd = Form("work-%08d%06d", td.GetDate(), td.GetTime());
  }
  gSystem->mkdir(wd);
  gSystem->CopyFile("ConfigFemtoAnalysis.C", wd + "/ConfigFemtoAnalysis.C");
  gSystem->cd(wd);

#if true
  //std::string file = "files-246001.txt";
  //std::string ffile = "files-AOD194.txt";
  std::string ffile = "../localfiles.txt";

  std::set<TString> input_files = load_files(ffile);

  TChain *input = new TChain("aodTree");

  for (auto &file : input_files) {
    std::cout << "Adding: " << file << "\n";
    input->Add(file);
  }

  mgr->InitAnalysis();
  mgr->PrintStatus();

  mgr->StartAnalysis("local", input, 350);
#else
  //wd = "AnalysisResult-12";
  setup_grid(mgr, wd);

  mgr->InitAnalysis();
  mgr->PrintStatus();

  mgr->StartAnalysis("grid");
#endif

  timer.Stop();
  timer.Print();

  TString outfile = wd + "/" + mgr->GetCommonFileName();
  std::cout << "Output written to " << outfile << "\n";
}

#if false

TChain* CreateChain(const char *xmlfile, const char *type="ESD");

const char *anatype = "AOD";

void myAnalysis()
{
// Analysis using AOD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// connect to AliEn and make the chain
   if (!TGrid::Connect("alien://")) return;
// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// load base root libraries
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");

   gSystem->Load("libMinuit");

// Load analysis framework libraries
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libOADB");
   gSystem->Load("libCORRFW");

// include path
   TString intPath = gInterpreter->GetIncludePath();
   TObjArray *listpaths = intPath.Tokenize(" ");
   TIter nextpath(listpaths);
   TObjString *pname;
   while ((pname=(TObjString*)nextpath())) {
      TString current = pname->GetName();
      if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
      gSystem->AddIncludePath(current);
   }
   if (listpaths) delete listpaths;
   gSystem->AddIncludePath("-I$ALICE_PHYSICS/include ");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());

// Add aditional AliRoot libraries

// analysis source to be compiled at runtime (if any)

// read the analysis manager from file
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("analysis.root");
   if (!mgr) return;
   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   TChain *chain = CreateChain("wn.xml", anatype);

   mgr->StartAnalysis("localfile", chain, 1234567890, 0);
   timer.Stop();
   timer.Print();
}

//________________________________________________________________________________
TChain* CreateChain(const char *xmlfile, const char *type)
{
// Create a chain using url's from xml file
   TString filename;
   Int_t run = 0;
   TString treename = type;
   treename.ToLower();
   treename += "Tree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   TAlienCollection *coll = dynamic_cast<TAlienCollection *>(TAlienCollection::Open(xmlfile));
   if (!coll) {
      ::Error("CreateChain", "Cannot create an AliEn collection from %s", xmlfile);
      return NULL;
   }
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TChain *chain = new TChain(treename);
   coll->Reset();
   while (coll->Next()) {
      filename = coll->GetTURL();
      if (mgr) {
         Int_t nrun = AliAnalysisManager::GetRunFromAlienPath(filename);
         if (nrun && nrun != run) {
            printf("### Run number detected from chain: %d\n", nrun);
            mgr->SetRunFromPath(nrun);
            run = nrun;
         }
      }
      chain->Add(filename);
   }
   if (!chain->GetNtrees()) {
      ::Error("CreateChain", "No tree found from collection %s", xmlfile);
      return NULL;
   }
   return chain;
}

#endif
