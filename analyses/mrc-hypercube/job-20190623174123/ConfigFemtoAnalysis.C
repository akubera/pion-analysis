///
/// \file mrc-hypercube/ConfigFemtoAnalysis.C
///

#if __cplusplus < 201103L
#error "This file requires use of the c++11 standard (ROOT6)"
#endif

#include <AliFemtoManager.h>
#include <AliFemtoEventReaderAlt.h>
#include <AliFemtoAnalysisPionPion.h>
#include <AliFemtoModelCorrFctnTrueQ6D.h>

using AFAPP = AliFemtoAnalysisPionPion;


void
AddEventReader(AliFemtoManager &mgr)
{
  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAlt();

  auto multest = AliFemtoEventReaderAOD::EstEventMult::kCentrality;
  const ULong_t filter_mask = BIT(7);

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

  auto *hyc_cf = new AliFemtoModelCorrFctnTrueQ6D::Builder()
                  .QoutRange(0.0, 0.1325)
                  .QsideRange(-0.1325, 0.1325)
                  .QlongRange(-0.1325, 0.1325)
                  .IntoCF();

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
  ccfg.event_CentralityMin = 0.0;
  ccfg.event_CentralityMax = 90.0;

  AFAPP::AnalysisParams acfg = AFAPP::DefaultConfig();
  acfg.num_events_to_mix = 6;
  acfg.enable_pair_monitors = false;
  acfg.is_mc_analysis = true;
  acfg.calc_automult(ccfg);

  acfg.pion_type_1 = AFAPP::kPiPlus;
  AddAnalysis("AnalysisMrc_pip", acfg, ccfg, *manager);

  acfg.pion_type_1 = AFAPP::kPiMinus;
  AddAnalysis("AnalysisMrc_pim", acfg, ccfg, *manager);

  return manager;
}
