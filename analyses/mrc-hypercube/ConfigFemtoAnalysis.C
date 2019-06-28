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

//static const ULong_t filter_mask = BIT(7);
static const ULong_t filter_mask = BIT(8);
// static const ULong_t filter_mask = BIT(5) | BIT(6);

void
AddEventReader(AliFemtoManager &mgr)
{
  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAlt();

  auto multest = AliFemtoEventReaderAOD::EstEventMult::kCentrality;

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
  analysis->AddStanardCutMonitors();

  auto *track_cut = new PrimaryPionCut();
  // analysis->SetFirstParticleCut(track_cut);
  // analysis->SetSecondParticleCut(track_cut);
  track_cut->ncls_its_min = 4;
  track_cut->pt_range = {0.14, 2.0};
  track_cut->eta_range = {-0.4, 0.4};
  track_cut->nsigma_pion = 1.0;
  track_cut->max_xy = 0.05;
  track_cut->max_z = 0.05;

  auto *ff = new AliFemtoModelManager();
  auto *wg = new AliFemtoModelWeightGeneratorBasic();
  ff->AcceptWeightGenerator(wg);

  auto *mrc_cf_1d = new AliFemtoModelCorrFctn("TrueQinv", 200, 0.0, 1.0);
  mrc_cf_1d->ConnectToManager(ff);

  auto *ktmrc1d = new AliFemtoKtBinnedCorrFunc("KT_MRC1D", mrc_cf_1d);


  // analysis->AddCorrFctn(mrc_cf_1d);

  auto *mrc_cf = new AliFemtoModelCorrFctnTrueQ6D("MRC",
                  // 22, 0.0, 0.1125,
                  // 45, -0.1125, 0.1125,
                  // 45, -0.1125, 0.1125);
                  26, 0.0, 0.1325,
                  53, -0.1325, 0.1325,
                  53, -0.1325, 0.1325);
  auto *ktmrc = new AliFemtoKtBinnedCorrFunc("KT_HYPERCUBE", mrc_cf);
#if true
  const std::vector<std::pair<float, float>> ktbins = {{0.2, 0.3}, {0.4, 0.5}, {0.6, 0.7}, {1.0, 1.2}};
#else
  const std::vector<std::pair<float, float>> ktbins = {{0.3, 0.4}, {0.5, 0.6}, {0.7, 0.8}, {0.8, 1.0}};
#endif

  ktmrc1d->AddKtRanges(ktbins);
  ktmrc->AddKtRanges(ktbins);

  analysis->AddCorrFctn(ktmrc1d);
  analysis->AddCorrFctn(ktmrc);
  m.AddAnalysis(analysis);
}

AliFemtoManager*
ConfigFemtoAnalysis()
{
  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AddEventReader(*manager);

  AFAPP::CutParams ccfg;
  ccfg.cuts_use_attrs = false;
  ccfg.mc_pion_only = true;
  ccfg.pion_1_rm_neg_lbl = true;
  ccfg.event_CentralityMin = 0.0;
  ccfg.event_CentralityMax = 90.0;

  AFAPP::AnalysisParams acfg = AFAPP::DefaultConfig();
  acfg.pion_type_2 = AFAPP::kNone;
  acfg.num_events_to_mix = 7;
  acfg.enable_pair_monitors = false;
  acfg.is_mc_analysis = true;
  // acfg.calc_automult(ccfg);

  acfg.pion_type_1 = AFAPP::kPiPlus;
  AddAnalysis("AnalysisMrc_pip", acfg, ccfg, *manager);

  // acfg.pion_type_1 = AFAPP::kPiMinus;
  // AddAnalysis("AnalysisMrc_pim", acfg, ccfg, *manager);

  return manager;
}
