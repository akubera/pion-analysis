///
/// \file ConfigFemtoAnalysis.C
///

#include <AliFemtoManager.h>
#include <AliFemtoEventReaderAlt.h>
#include <AliFemtoAnalysisPionPion.h>

#include <AliFemtoCorrFctn.h>
#include <AliFemtoCorrFctn3DLCMSPosQuad.h>

#include <AliFemtoAnalysisPionPionCuts.h>

#include <sstream>

using AFAPP = AliFemtoAnalysisPionPion;


struct PairPlotterCut : public AliFemtoPairCutPionPionAKDetaDphi {

  PairPlotterCut(const AliFemtoPairCutPionPionAKDetaDphi &orig)
    : AliFemtoPairCutPionPionAKDetaDphi(orig)
    {
    }

  bool Pass(const AliFemtoPair *pair) override
    {
      const AliFemtoTrack
        &track1 = *pair->Track1()->Track(),
        &track2 = *pair->Track2()->Track();

      if (!pwgfemto::PairCutAttrsBaseAK::Pass(*pair)) {
        return false;
      }

      if (PairCutTrackAttrDetaDphiStar::Pass(track1, track2)) {
          std::stringstream s1,
                            s2;

          for (int i=0; i<9; ++i) {
            const auto
              &p1 = track1.NominalTpcPoint(i),
              &p2 = track2.NominalTpcPoint(i);

            s1 << (i == 0 ? "" : ", ") << Form("[%g, %g, %g]", p1.x(), p1.y(), p1.z());
            s2 << (i == 0 ? "" : ", ") << Form("[%g, %g, %g]", p2.x(), p2.y(), p2.z());
          }

        std::cout << "[" << s1.str() << "]\n";
        std::cout << "[" << s2.str() << "]\n\n";

      }

      return true;
    }
};

/*
struct AliFoo : public AliFemtoCorrFctn {

  AliFoo()
    : AliFemtoCorrFctn()
    {
      std::cout << "AliFoo\n";
    }

  void AddMixedPair(AliFemtoPair *pair) override
    {
      AddPair(pair);
    }

  void AddRealPair(AliFemtoPair *pair) override
    {
      AddPair(pair);
    }

  void AddPair(AliFemtoPair *pair)
    {
      const double
        qo = pair->QOutCMS(),
        qs = pair->QSideCMS(),
        ql = pair->QLongCMS();

      const AliFemtoTrack
        &track1 = *pair->Track1()->Track(),
        &track2 = *pair->Track2()->Track();

      if (std::abs(ql) < 0.002) {
        if (std::abs(qo - 0.04) < 4e-3 && std::abs(qs - 0.00375) < 0.00375) {
          std::stringstream s1,
                            s2;

          for (int i=0; i<9; ++i) {
            const auto
              &p1 = track1.NominalTpcPoint(i),
              &p2 = track2.NominalTpcPoint(i);

            s1 << (i == 0 ? "" : ", ") << Form("[%g, %g, %g]", p1.x(), p1.y(), p1.z());
            s2 << (i == 0 ? "" : ", ") << Form("[%g, %g, %g]", p2.x(), p2.y(), p2.z());
          }

        std::cout << "[" << s1.str() << "]\n";
        std::cout << "[" << s2.str() << "]\n\n";
        }
      }
    }

  std::string Report() override
    {
      return "";
    }

  void Finish() override
    {
    }

  // void CreateOutputObjects() override
    // {}

  TList* GetOutputList() override
    {
      auto olist = new TList();
      return olist;
    }

  AliFoo* Clone() const override
    { return new AliFoo(); }

};
*/


void
AddEventReader(AliFemtoManager &mgr)
{
  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAlt();

  auto multest = AliFemtoEventReaderAOD::EstEventMult::kCentrality;

  rdr->SetFilterMask(BIT(7));
  rdr->SetEPVZERO(0);
  rdr->SetUseMultiplicity(multest);
  rdr->SetCentralityFlattening(false);
  rdr->SetReadV0(false);
  rdr->SetPrimaryVertexCorrectionTPCPoints(true);
  rdr->SetDCAglobalTrack(1);
  rdr->SetReadMC(false);
  rdr->SetReadFullMCData(false);
  mgr.SetEventReader(rdr);
}

void
AddPionAnalysis(AliFemtoManager &mgr)
{
  AFAPP::CutParams ccfg;
  // ccfg.cuts_use_attrs = true;
  // ccfg.mc_pion_only = true;
  // ccfg.pion_1_rm_neg_lbl = true;
  ccfg.pion_1_status = 16;
  // ccfg.pion_1_max_impact_xy = 0.023615054202417968;
  // ccfg.pion_1_max_impact_z = 0.03352419215870559;
  ccfg.pair_delta_eta_min = 0.0;
  ccfg.pair_delta_phi_min = 0.0;

  AFAPP::AnalysisParams acfg;
  acfg.is_mc_analysis = false;
  acfg.pion_type_1 = AFAPP::kPiPlus;
  acfg.pion_type_2 = AFAPP::kNone;

  auto *analysis = new AliFemtoAnalysisPionPion("pion_analysis_10_20_pip", acfg, ccfg);

  auto *pair_cut = new PairPlotterCut(*static_cast<AliFemtoPairCutPionPionAKDetaDphi*>(analysis->PairCut()));
  analysis->SetPairCut(pair_cut);

  // analysis->
  // analysis->AddCorrFctn(new AliFoo());
  mgr.AddAnalysis(analysis);
}

AliFemtoManager*
ConfigFemtoAnalysis()
{
  std::cout << "[ConfigFemtoAnalysis]\n";

  // Begin to build the manager and analyses
  auto *manager = new AliFemtoManager();

  AddEventReader(*manager);
  AddPionAnalysis(*manager);
  // exit(0);

  return manager;
}
