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

  PrimaryPionCut(const AliFemtoTrackCutPionPionIdealAK &cut)
    : AliFemtoTrackCutPionPionIdealAK(cut)
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

struct PseudoEsdTrackCut : public AliFemtoTrackCutPionPionIdealAK, public AliFemtoESDTrackCut {

  PseudoEsdTrackCut()
    : AliFemtoTrackCutPionPionIdealAK()
    , AliFemtoESDTrackCut()
    {}

  PseudoEsdTrackCut(AliFemtoConfigObject &cfg)
    : AliFemtoTrackCutPionPionIdealAK(cfg)
    , AliFemtoESDTrackCut()
    {}

  PseudoEsdTrackCut(const AliFemtoESDTrackCut &cut)
    : AliFemtoTrackCutPionPionIdealAK()
    , AliFemtoESDTrackCut(cut)
    {
      max_xy = fMaxImpactXY;
      max_z = fMaxImpactZ;
      pt_range = {fPt[0], fPt[1]};
      nsigma_pion = fNsigma;
      rchi2_tpc_max = fMaxTPCchiNdof;
      // rchi2_its_max = fMaxITSchiNdof;
      rchi2_its_max = 1.6;
    }

  virtual ~PseudoEsdTrackCut() {}

  virtual bool Pass(const AliFemtoTrack *track)
    {
      // const AliFemtoModelHiddenInfo *info = static_cast<AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
      // const Int_t origin = info->GetOrigin();
      // if (origin != 0) { // } || std::fabs(info->GetMass() - 0.139) > .001) {
      //   return false;
      // }

      if (!AliFemtoTrackCutPionPionAK::Pass(track)) {
        return false;
      }

      if (!Pass1(track)) {
        return false;
      }

      // if (!Pass2(track)) {
      //   return false;
      // }

      // if (!Pass3(track)) {
      //   return false;
      // }

      fNTracksPassed++;
      return true;
      // return AliFemtoESDTrackCut::Pass(track);
    }

  bool Pass1(const AliFemtoTrack *track)
    {
        if (fStatus && (track->Flags() & fStatus) != fStatus) {
          return false;
        }
        if (fRemoveKinks && (track->KinkIndex(0) || track->KinkIndex(1) || track->KinkIndex(2))) {
          return false;
        }
        if (fRemoveITSFake && track->ITSncls() < 0) {
          return false;
        }
        if (fminTPCclsF > track->TPCnclsF()) {
          return false;
        }
        if (fminTPCncls > track->TPCncls()) {
          return false;
        }
        if (fminITScls > track->ITSncls()) {
          return false;
        }

        if (fMaxImpactXY < TMath::Abs(track->ImpactD())) {
          return false;
        }
        if (fMinImpactXY > TMath::Abs(track->ImpactD())) {
          return false;
        }
        if (fMaxImpactZ < TMath::Abs(track->ImpactZ())) {
          return false;
        }
        if (fMaxSigmaToVertex < track->SigmaToVertex()) {
          return false;
        }

      /*
        if (track->ITSncls() > 0 && (track->ITSchi2() / track->ITSncls()) > fMaxITSchiNdof) {
          return false;
        }
        */

        if (track->TPCchi2perNDF() > fMaxTPCchiNdof) {
          return false;
        }

        // ITS cluster requirenments
        for (Int_t i = 0; i < 3; i++) {
          if (!CheckITSClusterRequirement(fCutClusterRequirementITS[i], track->HasPointOnITSLayer(i * 2), track->HasPointOnITSLayer(i*2+1))) {
            return false;
          }
        }

        if (fLabel) {
          if (track->Label() < 0) {
            fNTracksFailed++;
            return false;
          }
        }
        if (fCharge != 0 && (track->Charge() != fCharge)) {
          fNTracksFailed++;
          return false;
        }


        Bool_t tTPCPidIn = (track->Flags() & AliFemtoTrack::kTPCpid) > 0;
        Bool_t tITSPidIn = (track->Flags() & AliFemtoTrack::kITSpid) > 0;
        Bool_t tTOFPidIn = (track->Flags() & AliFemtoTrack::kTOFpid) > 0;

        const double momentum = track->P().Mag();

        if (fMinPforTOFpid > 0
            && fMinPforTOFpid < momentum && momentum < fMaxPforTOFpid
            && !tTOFPidIn) {
          fNTracksFailed++;
          return false;
        }

        if (fMinPforTPCpid > 0
            && fMinPforTPCpid < momentum && momentum < fMaxPforTPCpid
            && !tTPCPidIn) {
          fNTracksFailed++;
          return false;
        }

        if (fMinPforITSpid > 0
            && fMinPforITSpid < momentum && momentum < fMaxPforITSpid
            && !tITSPidIn) {
          fNTracksFailed++;
          return false;
        }

        float tEnergy = track->P().MassHypothesis(static_cast<AliFemtoESDTrackCut*>(this)->Mass());
        float tRapidity = 0;
        if (tEnergy-track->P().z() != 0 && (tEnergy + track->P().z()) / (tEnergy-track->P().z()) > 0)
          tRapidity = 0.5 * ::log((tEnergy + track->P().z())/(tEnergy-track->P().z()));
        float tPt = track->P().Perp();
        float tEta = track->P().PseudoRapidity();

        if (fMaxImpactXYPtOff < 999.0) {
          if ((fMaxImpactXYPtOff + fMaxImpactXYPtNrm*TMath::Power(tPt, fMaxImpactXYPtPow)) < TMath::Abs(track->ImpactD())) {
            fNTracksFailed++;
            return false;
          }
        }

        if ((tRapidity < fRapidity[0]) || (tRapidity > fRapidity[1])) {
          fNTracksFailed++;
          return false;
        }
        if ((tEta < fEta[0]) || (tEta > fEta[1])) {
          fNTracksFailed++;
          return false;
        }
        if ((tPt < fPt[0]) || (tPt > fPt[1])) {
          fNTracksFailed++;
          return false;
        }

      return true;
    }

  bool Pass2(const AliFemtoTrack *track)
    {
      if ((track->PidProbElectron() < fPidProbElectron[0]) || (track->PidProbElectron() > fPidProbElectron[1])) {
        fNTracksFailed++;
        return false;
      }
      if ((track->PidProbPion() < fPidProbPion[0]) || (track->PidProbPion() > fPidProbPion[1])) {
        fNTracksFailed++;
        return false;
      }
      if ((track->PidProbKaon() < fPidProbKaon[0]) || (track->PidProbKaon() > fPidProbKaon[1])) {
        fNTracksFailed++;
        return false;
      }
      if ((track->PidProbProton() < fPidProbProton[0]) || (track->PidProbProton() > fPidProbProton[1])) {
        fNTracksFailed++;
        return false;
      }
      if ((track->PidProbMuon() < fPidProbMuon[0]) || (track->PidProbMuon() > fPidProbMuon[1])) {
        fNTracksFailed++;
        return false;
      }
      return true;
    }

  bool Pass3(const AliFemtoTrack *track)
    {
      if (fElectronRejection && !IsElectron(track->NSigmaTPCE(),track->NSigmaTPCPi(),track->NSigmaTPCK(), track->NSigmaTPCP())) {
          return false;
      }

      if (!IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi())) {
        return false;
      }

      return true;
    }

  virtual const char* ClassName() const
    { return "PseudoEsdTrackCut"; }
};

// static const ULong_t filter_mask = BIT(5) | BIT(6);
static const ULong_t filter_mask = BIT(7);
// static const ULong_t filter_mask = BIT(8);

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


  // Use ESD Track cut
  // c.cuts_use_attrs = false;
  // auto *track_cut = analysis->BuildPionCut1(c);

  // Use Copy of ESD Track cut
  // c.cuts_use_attrs = false;
  // AliFemtoESDTrackCut *track_cut = new AliFemtoESDTrackCut(static_cast<const AliFemtoESDTrackCut&>(*analysis->BuildPionCut1(c)));

  // Use ESD Track cut with wide chi2-ITS
  c.cuts_use_attrs = false;
  auto *track_cut = static_cast<AliFemtoESDTrackCut *>(analysis->BuildPionCut1(c));
  track_cut->SetMaxITSChiNdof(1000);

  // Use Pseudo-ESD Track Cut
  // c.cuts_use_attrs = false;
  // AliFemtoESDTrackCut *track_cut = new PseudoEsdTrackCut(static_cast<const AliFemtoESDTrackCut&>(*analysis->BuildPionCut1(c)));

  // Use PrimaryPionCut
  // auto *track_cut = new PrimaryPionCut(*static_cast<AliFemtoTrackCutPionPionIdealAK*>(analysis->FirstParticleCut()));
  // track_cut->status = 16;


  analysis->SetFirstParticleCut(track_cut);
  analysis->SetSecondParticleCut(track_cut);
  // track_cut->ncls_its_min = 4;
  // track_cut->pt_range = {0.14, 2.0};
  // track_cut->eta_range = {-0.4, 0.4};
  // track_cut->nsigma_pion = 1.0;
  // track_cut->max_xy = 0.05;
  // track_cut->max_z = 0.05;

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
  ccfg.pion_1_rm_neg_lbl = false;
  ccfg.pion_1_status = 16;
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
