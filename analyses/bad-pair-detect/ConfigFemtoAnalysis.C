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


struct BadPairCut : public AliFemtoPairCutPionPionAKDetaDphi {

  BadPairCut(const AliFemtoPairCutPionPionAKDetaDphi &orig)
    : AliFemtoPairCutPionPionAKDetaDphi(orig)
    { }

  virtual bool Pass(const AliFemtoPair *pair)
    {
      // return AliFemtoPairCutPionPionAKDetaDphi::Pass(pair);
      const auto
        &track1 = *pair->Track1()->Track(),
        &track2 = *pair->Track2()->Track();

      auto is_selected_track = [] (const AliFemtoTrack &track)
        {
          auto *info = static_cast<AliFemtoModelHiddenInfo*>(track.GetHiddenInfo());
          if (info->GetPDGPid() == -211) {
            return true;
          }
          // double p = track.P().Mag();
          // if (std::fabs(p - 0.14589) < 6e-5) {
          //   return true;
          // }
          return false;
        };

      if (is_selected_track(track1) || is_selected_track(track2)) {
        return AliFemtoPairCutPionPionAKDetaDphi::Pass(pair);
      }

      return false;
    }

  virtual const char* ClassName() const
    {
      return "BadPairCut";
    }
};


struct BadPairFinder : public AliFemtoModelCorrFctn {
  BadPairFinder()
    : AliFemtoModelCorrFctn()
    {
  delete fNumeratorTrue;
  fNumeratorTrue = nullptr;
  delete fNumeratorFake;
  fNumeratorFake = nullptr;
  delete fDenominator;
  fDenominator = nullptr;
  delete fNumeratorTrueIdeal;
  fNumeratorTrueIdeal = nullptr;
  delete fNumeratorFakeIdeal;
  fNumeratorFakeIdeal = nullptr;
  delete fDenominatorIdeal;
  fDenominatorIdeal = nullptr;

  delete fQgenQrec;
  fQgenQrec = nullptr;
    }

  virtual ~BadPairFinder()
    { }

  void AddPair(const AliFemtoPair &pair)
    {
      const auto
        &track1 = *pair.Track1(),
        &track2 = *pair.Track2();

      const auto
        *info1 = (AliFemtoModelHiddenInfo*)track1.GetHiddenInfo(),
        *info2 = (AliFemtoModelHiddenInfo*)track2.GetHiddenInfo();

      const auto
        &mp1 = *info1->GetTrueMomentum(),
        &mp2 = *info2->GetTrueMomentum(),

        *op1 = track1.Track()->GetTrueMomentum(),
        *op2 = track2.Track()->GetTrueMomentum();

      const AliFemtoLorentzVector
        tp1(mp1, mp1.MassHypothesis(info1->GetMass())),
        tp2(mp2, mp2.MassHypothesis(info2->GetMass())),

        &p1 = track1.FourMomentum(),
        &p2 = track2.FourMomentum();

      const double
        qinv = - (p1 - p2).m(),
        true_qinv = - (tp1 - tp2).m();

      const auto
        &mom1 = track1.Track()->P(),
        &mom2 = track2.Track()->P();

      if (true_qinv > 0.4 and qinv <= (true_qinv - 0.3)) {
        std::cout << "Found bad pair!"
                  << "\n"
                  << Form(" q: (%0.5f, %0.5f)", true_qinv, qinv)
                  // << Form(" truep: (%0.5f, %0.5f)", mp1.Mag(), mp2.Mag())
                  << "\n"
                  << "    p: " << Form("<%g,%g,%g: %g>", mom1.x(), mom1.y(), mom1.z(), mom1.Mag())
                    << " "  << Form("<%g,%g,%g: %g>", mom2.x(), mom2.y(), mom2.z(), mom2.Mag())
                  << "\n"
                  << " true: " << Form("<%g,%g,%g: %g>", mp1.x(), mp1.y(), mp1.z(), mp1.Mag())
                    << " "  << Form("<%g,%g,%g: %g>", mp2.x(), mp2.y(), mp2.z(), mp2.Mag())
                  // << " truep: " << mp1 << "," << mp2
                  << " otruep: (" << (op1 ? Form("%0.5f", op1->Mag()) : "NULL") << ", " << (op2 ? Form("%0.5f", op2->Mag()) : "NULL") << ")"
                  // << Form(" otruep: (%0.5f, %0.5f)", track1.Track()->GetTrueMomentum()->Mag(), track2.Track()->GetTrueMomentum()->Mag())
                  // << Form(" p: (%0.5f, %0.5f)", track1.Track()->P().Mag(), track2.Track()->P().Mag())
                  << "\n"
                  << "  PDG: " << info1->GetPDGPid() << ", " << info2->GetPDGPid()
                  << " Mother: " << info1->GetMotherPdgCode() << ", " << info2->GetMotherPdgCode()
                  << " Origin: " << info1->GetOrigin() << ", " << info2->GetOrigin()
                  << " TpcChi2:" << track1.Track()->TPCchi2perNDF() << ", " << track2.Track()->TPCchi2perNDF()
                  << " ItsChi2:" << track1.Track()->ITSchi2perNDF() << ", " << track2.Track()->ITSchi2perNDF() << "\n"
                  << "\n";
      }
    }

  void AddRealPair(AliFemtoPair *pair) override
    {
      AddPair(*pair);
    }

  void AddMixedPair(AliFemtoPair *pair) override
    {
      AddPair(*pair);
    }

  TList* GetOutputList() override
    {
      return new TList();
    }
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

  auto *pair_cut = new BadPairCut(*static_cast<AliFemtoPairCutPionPionAKDetaDphi*>(analysis->PairCut()));
  analysis->SetPairCut(pair_cut);

  auto *ff = new AliFemtoModelManager();
  auto *wg = new AliFemtoModelWeightGeneratorBasic();
  ff->AcceptWeightGenerator(wg);

  auto *mrc_cf_1d = new AliFemtoModelCorrFctn("TrueQinv", 350, 0.0, 1.2);
  mrc_cf_1d->ConnectToManager(ff);
  analysis->AddCorrFctn(mrc_cf_1d);

  auto *badpair_cf = new BadPairFinder();
  analysis->AddCorrFctn(badpair_cf);

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
  ccfg.pion_1_PtMin = 0.10;
  ccfg.pion_1_min_its_ncls = 3;
  ccfg.pion_1_max_its_chi_ndof = 2.50;
  ccfg.pion_1_max_tpc_chi_ndof = 2.60;
  ccfg.pion_1_min_tpc_chi_ndof = 0.30;

  AFAPP::AnalysisParams acfg = AFAPP::DefaultConfig();
  acfg.pion_type_2 = AFAPP::kNone;
  acfg.num_events_to_mix = 12;
  acfg.enable_pair_monitors = false;
  acfg.is_mc_analysis = true;
  // acfg.calc_automult(ccfg);

  acfg.pion_type_1 = AFAPP::kPiPlus;
  AddAnalysis("AnalysisMrc_pip", acfg, ccfg, *manager);

  // acfg.pion_type_1 = AFAPP::kPiMinus;
  // AddAnalysis("AnalysisMrc_pim", acfg, ccfg, *manager);

  return manager;
}
