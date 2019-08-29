///
/// \file mrc-hypercube/ConfigFemtoAnalysis.C
///

#if __cplusplus < 201103L
#error "This file requires use of the c++11 standard (ROOT6)"
#endif

#include <AliFemtoManager.h>
#include <AliFemtoEventReaderAlt.h>
#include <AliFemtoAnalysisPionPion.h>
#include <AliFemtoModelCorrFctnQinv.h>
#include <AliFemtoModelCorrFctnTrueQ6D.h>
#include <AliFemtoKtBinnedCorrFunc.h>
#include <AliFemtoAnalysisPionPionCuts.h>
#include <AliFemtoModelManager.h>
#include <AliFemtoModelWeightGeneratorBasic.h>
#include <AliFemtoESDTrackCut.h>

using AFAPP = AliFemtoAnalysisPionPion;


struct PrimaryPionCut : public AliFemtoTrackCutPionPionIdealAK {

  PrimaryPionCut()
    : AliFemtoTrackCutPionPionIdealAK()
    {}

  PrimaryPionCut(AliFemtoConfigObject &cfg)
    : AliFemtoTrackCutPionPionIdealAK(cfg)
    {}

  virtual ~PrimaryPionCut() {}

  virtual bool Pass(const AliFemtoTrack *track)
    {
      const AliFemtoModelHiddenInfo *info = static_cast<AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
      const Int_t origin = info->GetOrigin();
      if (origin != 0 || std::fabs(info->GetMass() - 0.139) > .001) {
        return false;
      }
      return AliFemtoTrackCutPionPionIdealAK::Pass(track);
    }

  virtual const char* ClassName() const
    { return "PrimaryPionCut"; }
};


static const ULong_t filter_mask = BIT(7);
// static const ULong_t filter_mask = BIT(8);
// static const ULong_t filter_mask = BIT(5) | BIT(6);

void
AddEventReader(AliFemtoManager &mgr)
{
  auto *rdr = new AliFemtoEventReaderAlt();

  auto multest = AliFemtoEventReaderAOD::EstEventMult::kCentrality;

  rdr->SetEnhanceSmearing(0.1);

  rdr->SetFilterMask(filter_mask);
  rdr->SetEPVZERO(0);
  rdr->SetUseMultiplicity(multest);
  rdr->SetCentralityFlattening(false);
  rdr->SetReadV0(false);
  rdr->SetPrimaryVertexCorrectionTPCPoints(true);
  rdr->SetDCAglobalTrack(1);
  rdr->SetReadMC(true);
  rdr->SetReadFullMCData(true);
  mgr.SetEventReader(rdr);
}


void
AddAnalysis(TString name, AFAPP::AnalysisParams a, AFAPP::CutParams c, AliFemtoManager &m)
{
  AliFemtoAnalysisPionPion *analysis = new AliFemtoAnalysisPionPion(name, a, c);
  analysis->SetTrackFilter(filter_mask);
/*
  c.cuts_use_attrs = false;
  auto *track_cut = static_cast<AliFemtoESDTrackCut*>(analysis->BuildPionCut1(c));
  c.cuts_use_attrs = true;

  //auto *track_cut = new PrimaryPionCut();
  // auto *track_cut = new AliFemtoTrackCutPionPionIdealAK();
  analysis->SetFirstParticleCut(track_cut);
  analysis->SetSecondParticleCut(track_cut);
*/

#if true
  const std::vector<std::pair<float, float>> ktbins = {{0.2, 0.3}, {0.4, 0.5}, {0.6, 0.7}};
#else
  const std::vector<std::pair<float, float>> ktbins = {{0.3, 0.4}, {0.5, 0.6}, {0.7, 0.8}, {0.8, 1.0}, {1.0, 1.2}};
#endif

  auto *ff = new AliFemtoModelManager();
  auto *wg = new AliFemtoModelWeightGeneratorBasic();
  ff->AcceptWeightGenerator(wg);

  auto *mrc_cf_1d = new AliFemtoModelCorrFctnQinv("TrueQinv", 200, 0.0, 1.0);
  mrc_cf_1d->ConnectToManager(ff);

  auto *ktmrc1d = new AliFemtoKtBinnedCorrFunc("KT_MRC1D", mrc_cf_1d);
  ktmrc1d->AddKtRanges(ktbins);

  analysis->AddCorrFctn(ktmrc1d);

  // bins match CF_PbPb#6986
  int out_nbins = 38;
  double out_max = out_nbins * (0.1525000035762787 / 30);
  auto *mrc_cf = new AliFemtoModelCorrFctnTrueQ6D("MRC",
                  out_nbins, 0.0, out_max,
                  61, -0.1525000035762787, 0.1525000035762787,
                  61, -0.1525000035762787, 0.1525000035762787);

  auto *ktmrc = new AliFemtoKtBinnedCorrFunc("KT_HYPERCUBE", mrc_cf);

  ktmrc->AddKtRanges(ktbins);

  analysis->AddCorrFctn(ktmrc);

  analysis->AddStanardCutMonitors();

  m.AddAnalysis(analysis);
}


AliFemtoManager*
ConfigFemtoAnalysis()
{
  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AddEventReader(*manager);

  AFAPP::CutParams ccfg;
  ccfg.cuts_use_attrs = true;
  ccfg.mc_pion_only = true;
  ccfg.pion_1_rm_neg_lbl = true;
  ccfg.pion_1_status = 16;
  ccfg.pion_1_max_impact_xy = 0.03;
  ccfg.pion_1_max_impact_z = 0.04;
  ccfg.pion_1_sigma = 3.0;
  ccfg.pion_1_min_its_ncls = 3;
  ccfg.pion_1_max_its_chi_ndof = 2.5;
  ccfg.pion_1_max_tpc_chi_ndof = 1.6;
  ccfg.pion_1_min_tpc_chi_ndof = 0.33;
  ccfg.event_CentralityMin = 0.0;
  ccfg.event_CentralityMax = 90.0;

  AFAPP::AnalysisParams acfg = AFAPP::DefaultConfig();
  acfg.pion_type_2 = AFAPP::kNone;
  acfg.num_events_to_mix = 8;
  acfg.enable_pair_monitors = false;
  acfg.is_mc_analysis = true;
  // acfg.calc_automult(ccfg);

  acfg.pion_type_1 = AFAPP::kPiPlus;
  AddAnalysis("AnalysisMrc_pip", acfg, ccfg, *manager);

  acfg.pion_type_1 = AFAPP::kPiMinus;
  AddAnalysis("AnalysisMrc_pim", acfg, ccfg, *manager);

  return manager;
}
