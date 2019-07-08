///
/// \file mrc-1d/ConfigFemtoAnalysis.C
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

#include <TH2D.h>


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


class MrcMatrixCorrFctn : public AliFemtoCorrFctn {
protected:
  TH2D *fMatrixQinv;
  TH2D *fMatrixQlcms;

public:

  MrcMatrixCorrFctn()
    : AliFemtoCorrFctn()
    , fMatrixQinv(new TH2D("QgenQrec", "MRC Matrix; q_{inv,gen}; q_{inv,rec};", 400, 0.0, 1.0, 400, 0.0, 1.0))
    , fMatrixQlcms(new TH2D("QgenQreclcms", "MRC Matrix; q_{lcms,gen}; q_{lcms,rec};", 400, 0.0, 1.0, 400, 0.0, 1.0))
    {}

  MrcMatrixCorrFctn(const MrcMatrixCorrFctn &orig)
    : AliFemtoCorrFctn(orig)
    , fMatrixQinv(new TH2D(*orig.fMatrixQinv))
    , fMatrixQlcms(new TH2D(*orig.fMatrixQlcms))
    {}

  MrcMatrixCorrFctn& operator=(const MrcMatrixCorrFctn &rhs)
    {
      AliFemtoCorrFctn::operator=(rhs);
      *fMatrixQinv = *rhs.fMatrixQinv;
      *fMatrixQlcms = *rhs.fMatrixQlcms;
      return *this;
    }

  virtual MrcMatrixCorrFctn* Clone() const
    {
      return new MrcMatrixCorrFctn(*this);
    }

  static double CalcQinv(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
    {
      return AliFemtoModelCorrFctnQinv::CalcQinv(p1, p2);
    }

  static double CalcQlcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
    {
      const AliFemtoLorentzVector d = p1 - p2;

      const double
        beta = (p1.z() + p2.z()) / (p1.e() + p1.e()),

        x = (d.z() - beta * d.e()),

        gamma_sqrd = 1.0 / (1.0 - beta * beta),
        // gamma_sqrd = 1.0 / ((1.0 - beta) * (1.0 + beta)),

        qlong_sqrd = gamma_sqrd * x * x,

        qlcms = std::sqrt(d.Perp2() + qlong_sqrd);

      return qlcms;
    }

  virtual void AddRealPair(AliFemtoPair *)
    {}

  virtual void AddMixedPair(AliFemtoPair *pair)
    {
      const AliFemtoModelHiddenInfo
        *info1 = static_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->HiddenInfo()),
        *info2 = static_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->HiddenInfo());

      const AliFemtoLorentzVector
        &p1(pair->Track1()->FourMomentum()),
        &p2(pair->Track2()->FourMomentum()),
        ideal_p1(info1->GetTrueMomentum()->MassHypothesis(info1->GetMass()), *info1->GetTrueMomentum()),
        ideal_p2(info2->GetTrueMomentum()->MassHypothesis(info2->GetMass()), *info2->GetTrueMomentum());

      const double
        qinv_rec = CalcQinv(p1, p2),
        qinv_gen = CalcQinv(ideal_p1, ideal_p2),
        qlcms_rec = CalcQlcms(p1, p2),
        qlcms_gen = CalcQlcms(ideal_p1, ideal_p2);

      fMatrixQinv->Fill(qinv_gen, qinv_rec);
      fMatrixQlcms->Fill(qlcms_gen, qlcms_rec);
    }

  virtual void Finish()
    {}

  virtual AliFemtoString Report()
    {
      return "";
    }

  virtual TList* GetOutputList()
    {
      TList *list = new TList();
      list->Add(fMatrixQinv);
      list->Add(fMatrixQlcms);
      return list;
    }

};

static const ULong_t filter_mask = BIT(7);
// static const ULong_t filter_mask = BIT(8);
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

  // c.cuts_use_attrs = false;
  // auto *track_cut = static_cast<AliFemtoESDTrackCut*>(analysis->BuildPionCut1(c));
  // c.cuts_use_attrs = true;

  auto *track_cut = new PrimaryPionCut();
  // // auto *track_cut = new AliFemtoTrackCutPionPionIdealAK();
  analysis->SetFirstParticleCut(track_cut);
  analysis->SetSecondParticleCut(track_cut);

  // auto *ff = new AliFemtoModelManager();
  // auto *wg = new AliFemtoModelWeightGeneratorBasic();
  // ff->AcceptWeightGenerator(wg);

  auto *mrc_cf_1d = new MrcMatrixCorrFctn();

  auto *ktmrc1d = new AliFemtoKtBinnedCorrFunc("KT_MRC1D", mrc_cf_1d);

  // analysis->AddCorrFctn(mrc_cf_1d);

  // auto *mrc_cf = new AliFemtoModelCorrFctnTrueQ6D("MRC",
  //                 // 22, 0.0, 0.1125,
  //                 // 45, -0.1125, 0.1125,
  //                 // 45, -0.1125, 0.1125);
  //                 30, 0.0, 0.1525,
  //                 61, -0.1525, 0.1525,
  //                 61, -0.1525, 0.1525);

//   auto *ktmrc = new AliFemtoKtBinnedCorrFunc("KT_HYPERCUBE", mrc_cf);
// #if true
//   const std::vector<std::pair<float, float>> ktbins = {{0.2, 0.3}, {0.4, 0.5}, {0.6, 0.7}};
// #else
//   const std::vector<std::pair<float, float>> ktbins = {{0.3, 0.4}, {0.5, 0.6}, {0.7, 0.8}, {0.8, 1.0}, {1.0, 1.2}};
// #endif
  const std::vector<std::pair<float, float>> ktbins = {
    {0.2, 0.3}, {0.4, 0.5}, {0.6, 0.7},
    {0.3, 0.4}, {0.5, 0.6}, {0.7, 0.8}, {0.8, 1.0}, {1.0, 1.2}
  };

  ktmrc1d->AddKtRanges(ktbins);
  // ktmrc->AddKtRanges(ktbins);

  analysis->AddCorrFctn(ktmrc1d);
  // analysis->AddCorrFctn(ktmrc);

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
  ccfg.pion_1_max_impact_xy = 0.04;
  ccfg.pion_1_max_impact_z = 0.05;
  ccfg.pion_1_sigma = 2.5;
  ccfg.pion_1_min_its_ncls = 4;
  ccfg.pion_1_max_its_chi_ndof = 1.6;
  ccfg.pion_1_max_tpc_chi_ndof = 1.6;
  ccfg.pion_1_min_tpc_chi_ndof = 0.33;
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

  acfg.pion_type_1 = AFAPP::kPiMinus;
  AddAnalysis("AnalysisMrc_pim", acfg, ccfg, *manager);

  return manager;
}
