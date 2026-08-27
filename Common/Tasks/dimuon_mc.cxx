// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file dimuon_mc.cxx
/// \brief Pair-pT and rapidity resolution for reconstructed omega/phi -> mu+mu-.

#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{
using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::MostProbableEMEventIdsInMC>;
using MyMCMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds, aod::EMPrimaryMuonMCLabels, aod::EMMFTMCLabels>;
using FilteredMyMCMuons = soa::Filtered<MyMCMuons>;

struct TracksFromMother {
  std::vector<int64_t> positiveTrackIds;
  std::vector<int64_t> negativeTrackIds;
};

struct TracksFromEvent {
  std::vector<int64_t> positiveTrackIds;
  std::vector<int64_t> negativeTrackIds;
};
} // namespace

struct DimuonMcResolution {
  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "If non-negative, select this MC sub-generator ID"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "Centrality estimator: FT0M=0, FT0A=1, FT0C=2"};
  Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "Require both muons to belong to the MC event associated with the reconstructed event"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "Minimum reconstructed vertex z (cm)"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "Maximum reconstructed vertex z (cm)"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "Require sel8"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "Require FT0 coincidence"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "Require no time-frame border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "Require no ITS readout-frame border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "Require no same-bunch pileup"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "Require an ITS-TPC vertex contributor"};
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "Require a TOF-matched vertex contributor"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "Require compatible FT0 and PV z vertices"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "Minimum track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "Maximum track occupancy (exclusive)"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2.f, "Minimum FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000.f, "Maximum FT0C occupancy (exclusive)"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "Require no collision in the standard time range"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "Require no collision in the strict time range"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "Require no collision in the standard ITS ROF"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "Require no collision in the strict ITS ROF"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "Require no high-multiplicity collision in the previous ITS ROF"};
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "Require a good ITS layer 3"};
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "Require good ITS layers 0-3"};
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "Require all ITS layers to be good"};
    Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "Require the configured RCT detector-quality flag"};
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "RCT detector-quality label"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC in the RCT check"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "Treat limited detector acceptance as bad in the RCT check"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "Minimum centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "Maximum centrality"};
    Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "Minimum number of PV contributors"};
    Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "Maximum number of PV contributors (exclusive)"};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.f, "Minimum dimuon mass (GeV/c^2)"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1.e+10f, "Maximum dimuon mass (GeV/c^2)"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.f, "Minimum dimuon pT (GeV/c)"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1.e+10f, "Maximum dimuon pT (GeV/c)"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -4.f, "Minimum dimuon rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", -2.5f, "Maximum dimuon rapidity"};
    Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.f, "Minimum pair DCAxy"};
    Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1.e+10f, "Maximum pair DCAxy"};
    Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "Apply the pair delta-eta/delta-phi cut"};
    Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02f, "Minimum pair delta eta"};
    Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.02f, "Minimum pair delta phi"};
    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "Muon track type: 0 global, 3 MCH-MID standalone"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2f, "Minimum track pT (GeV/c)"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1.e+10f, "Maximum track pT (GeV/c)"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.f, "Minimum track eta"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5f, "Maximum track eta"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "Minimum track phi (rad)"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3f, "Maximum track phi (rad)"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "Minimum number of MFT clusters"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "Minimum number of MCH clusters"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1.e+6f, "Maximum track chi2/ndf"};
    Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1.e+6f, "Maximum MFT track chi2/ndf"};
    Configurable<float> cfg_border_pt_for_chi2mchmft{"cfg_border_pt_for_chi2mchmft", 0.f, "pT boundary for the MFT-MCH matching-chi2 cut (GeV/c)"};
    Configurable<float> cfg_max_matching_chi2_mftmch_lowPt{"cfg_max_matching_chi2_mftmch_lowPt", 8.f, "Maximum MFT-MCH matching chi2 below the pT boundary"};
    Configurable<float> cfg_max_matching_chi2_mftmch_highPt{"cfg_max_matching_chi2_mftmch_highPt", 40.f, "Maximum MFT-MCH matching chi2 above the pT boundary"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1.e+10f, "Maximum MCH-MID matching chi2"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.e+10f, "Maximum track DCAxy (cm)"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6f, "Minimum absorber-end radius (cm)"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5f, "Maximum absorber-end radius (cm)"};
    Configurable<float> cfg_max_diff_chi2_mftmch{"cfg_max_diff_chi2_mftmch", -1.f, "Maximum difference between best and second-best MFT-MCH matching chi2"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Enable track-to-collision association"};
    Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1.e+10f, "Maximum relative pT difference from the matched MCH-MID track"};
    Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1.e+10f, "Maximum eta difference from the matched MCH-MID track"};
    Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1.e+10f, "Maximum phi difference from the matched MCH-MID track"};
    Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "Require hits in configured MFT disks"};
    Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "Required MFT disks (logical OR)"};
    Configurable<bool> acceptOnlyCorrectMatch{"acceptOnlyCorrectMatch", false, "Accept only correctly matched MFT-MCH tracks"};
    Configurable<bool> acceptOnlyWrongMatch{"acceptOnlyWrongMatch", false, "Accept only incorrectly matched MFT-MCH tracks"};
    Configurable<bool> acceptOnlyXORMatching{"acceptOnlyXORMatching", false, "Accept global-muon pairs with exactly one correct MFT-MCH match"};
  } dimuoncuts;

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.1f, "Minimum generated-muon pT (GeV/c)"};
    Configurable<float> max_mcPt{"max_mcPt", 1.e+10f, "Maximum generated-muon pT (GeV/c)"};
    Configurable<float> min_mcEta{"min_mcEta", -0.8f, "Minimum generated-muon eta"};
    Configurable<float> max_mcEta{"max_mcEta", +0.8f, "Maximum generated-muon eta"};
  } mctrackcuts;

  ConfigurableAxis axisPtGen{"axisPtGen", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Generated meson pT binning"};
  ConfigurableAxis axisYGen{"axisYGen", {75, -4.f, -2.5f}, "Generated meson rapidity binning"};
  ConfigurableAxis axisRelPtResolution{"axisRelPtResolution", {200, -1.f, +1.f}, "Relative dimuon pT resolution binning"};
  ConfigurableAxis axisYResolution{"axisYResolution", {200, -0.5f, +0.5f}, "Dimuon rapidity resolution binning"};
  ConfigurableAxis axisMass{"axisMass", {240, 0.f, 1.2f}, "Reconstructed dimuon mass binning"};

  expressions::Filter trackFilter = aod::fwdtrack::trackType == dimuoncuts.cfg_track_type;
  expressions::Filter ttcaFilter = ifnode(dimuoncuts.enableTTCA.node(), true, aod::emprimarymuon::isAssociatedToMPC == true);
  expressions::Filter collisionFilterCentrality = (eventcuts.cfgCentMin < aod::cent::centFT0M && aod::cent::centFT0M < eventcuts.cfgCentMax) ||
                                                  (eventcuts.cfgCentMin < aod::cent::centFT0A && aod::cent::centFT0A < eventcuts.cfgCentMax) ||
                                                  (eventcuts.cfgCentMin < aod::cent::centFT0C && aod::cent::centFT0C < eventcuts.cfgCentMax);
  expressions::Filter collisionFilterNumContrib = eventcuts.cfgNumContribMin <= aod::collision::numContrib && aod::collision::numContrib < eventcuts.cfgNumContribMax;
  expressions::Filter collisionFilterTrackOccupancy = eventcuts.cfgTrackOccupancyMin <= aod::evsel::trackOccupancyInTimeRange && aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  expressions::Filter collisionFilterFT0COccupancy = eventcuts.cfgFT0COccupancyMin <= aod::evsel::ft0cOccupancyInTimeRange && aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;

  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  EMEventCut eventCut;
  DimuonCut dimuonCut;
  aod::rctsel::RCTFlagsChecker rctChecker;
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void addSpeciesHistograms(const std::string& species)
  {
    const std::string prefix = species + "/";
    registry.add(prefix + "hPtResolution", "Relative dimuon transverse-momentum resolution;(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen};pairs", HistType::kTH1F, {axisRelPtResolution});
    registry.add(prefix + "hYResolution", "Dimuon rapidity resolution;#it{y}^{gen} - #it{y}^{rec};pairs", HistType::kTH1F, {axisYResolution});
    registry.add(prefix + "hPtGenVsResolution", "Relative dimuon transverse-momentum resolution versus generated meson transverse momentum;#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}", HistType::kTH2F, {axisPtGen, axisRelPtResolution});
    registry.add(prefix + "hYGenVsResolution", "Dimuon rapidity resolution versus generated meson rapidity;#it{y}^{gen};#it{y}^{gen} - #it{y}^{rec}", HistType::kTH2F, {axisYGen, axisYResolution});
    registry.add(prefix + "hPtResponse", "Dimuon transverse-momentum response;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", HistType::kTH2F, {axisPtGen, axisPtGen});
    registry.add(prefix + "hYResponse", "Dimuon rapidity response;#it{y}^{gen};#it{y}^{rec}", HistType::kTH2F, {axisYGen, axisYGen});
    registry.add(prefix + "hMassVsPtGen", "Reconstructed dimuon mass versus generated meson transverse momentum;#it{p}_{T}^{gen} (GeV/#it{c});#it{m}_{#mu#mu}^{rec} (GeV/#it{c}^{2})", HistType::kTH2F, {axisPtGen, axisMass});
  }

  void defineEventCut()
  {
    eventCut = EMEventCut("eventCut", "eventCut");
    eventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    eventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    eventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    eventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    eventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    eventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    eventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    eventCut.SetRequireVertexTOFmatched(eventcuts.cfgRequireVertexTOFmatched);
    eventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    eventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    eventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    eventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    eventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    eventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
    eventCut.SetRequireGoodITSLayer3(eventcuts.cfgRequireGoodITSLayer3);
    eventCut.SetRequireGoodITSLayer0123(eventcuts.cfgRequireGoodITSLayer0123);
    eventCut.SetRequireGoodITSLayersAll(eventcuts.cfgRequireGoodITSLayersAll);
  }

  void defineDimuonCut()
  {
    dimuonCut = DimuonCut("dimuonCut", "dimuonCut");
    dimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    dimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    dimuonCut.SetPairYRange(dimuoncuts.cfg_min_pair_y, dimuoncuts.cfg_max_pair_y);
    dimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy);
    dimuonCut.SetMindEtadPhi(dimuoncuts.cfg_apply_detadphi, dimuoncuts.cfg_min_deta, dimuoncuts.cfg_min_dphi);
    dimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    dimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, dimuoncuts.cfg_max_pt_track);
    dimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    dimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    dimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    dimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 20);
    dimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    dimuonCut.SetChi2MFT(0.f, dimuoncuts.cfg_max_chi2mft);
    dimuonCut.SetMaxDiffMatchingChi2MCHMFT(dimuoncuts.cfg_max_diff_chi2_mftmch);
    dimuonCut.SetMaxMatchingChi2MCHMFTPtDep([&](float pt) { return pt < dimuoncuts.cfg_border_pt_for_chi2mchmft ? dimuoncuts.cfg_max_matching_chi2_mftmch_lowPt : dimuoncuts.cfg_max_matching_chi2_mftmch_highPt; });
    dimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    dimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    dimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    dimuonCut.SetMaxPDCARabsDep([](float rabs) { return rabs < 26.5f ? 594.f : 324.f; });
    dimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID);
    dimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
    dimuonCut.EnableTTCA(dimuoncuts.enableTTCA);
  }

  void init(InitContext const&)
  {
    if (cfgCentEstimator.value < 0 || cfgCentEstimator.value > 2) {
      LOGF(fatal, "cfgCentEstimator must be 0 (FT0M), 1 (FT0A), or 2 (FT0C)");
    }
    if (dimuoncuts.acceptOnlyCorrectMatch.value && dimuoncuts.acceptOnlyWrongMatch.value) {
      LOGF(fatal, "acceptOnlyCorrectMatch and acceptOnlyWrongMatch cannot both be enabled");
    }
    defineEventCut();
    defineDimuonCut();
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);
    addSpeciesHistograms("Omega");
    addSpeciesHistograms("Phi");
  }

  template <typename TTrack>
  bool isGeneratedMuonAccepted(TTrack const& track)
  {
    return mctrackcuts.min_mcPt.value < track.pt() && track.pt() < mctrackcuts.max_mcPt.value &&
           mctrackcuts.min_mcEta.value < track.eta() && track.eta() < mctrackcuts.max_mcEta.value;
  }

  template <typename TTrack1, typename TTrack2>
  bool isRequestedMatchingCategory(TTrack1 const& first, TTrack2 const& second)
  {
    const bool firstCorrect = first.emmcparticleId() == first.emmftmcparticleId();
    const bool secondCorrect = second.emmcparticleId() == second.emmftmcparticleId();
    const auto globalMuonType = static_cast<uint8_t>(aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack);

    if (dimuoncuts.acceptOnlyCorrectMatch.value &&
        ((first.trackType() == globalMuonType && !firstCorrect) || (second.trackType() == globalMuonType && !secondCorrect))) {
      return false;
    }
    if (dimuoncuts.acceptOnlyWrongMatch.value &&
        ((first.trackType() == globalMuonType && firstCorrect) || (second.trackType() == globalMuonType && secondCorrect))) {
      return false;
    }
    if (dimuoncuts.acceptOnlyXORMatching.value && first.trackType() == globalMuonType && second.trackType() == globalMuonType && !(firstCorrect ^ secondCorrect)) {
      return false;
    }
    return true;
  }

  template <int MotherPdg>
  void fillResolutionHistograms(float ptGen, float yGen, float ptRec, float yRec, float massRec)
  {
    const float yResolution = yGen - yRec;
    if constexpr (MotherPdg == 223) {
      registry.fill(HIST("Omega/hYResolution"), yResolution);
      registry.fill(HIST("Omega/hYGenVsResolution"), yGen, yResolution);
      registry.fill(HIST("Omega/hYResponse"), yGen, yRec);
      registry.fill(HIST("Omega/hMassVsPtGen"), ptGen, massRec);
      if (ptGen > 0.f) {
        const float ptResolution = (ptGen - ptRec) / ptGen;
        registry.fill(HIST("Omega/hPtResolution"), ptResolution);
        registry.fill(HIST("Omega/hPtGenVsResolution"), ptGen, ptResolution);
        registry.fill(HIST("Omega/hPtResponse"), ptGen, ptRec);
      }
    } else if constexpr (MotherPdg == 333) {
      registry.fill(HIST("Phi/hYResolution"), yResolution);
      registry.fill(HIST("Phi/hYGenVsResolution"), yGen, yResolution);
      registry.fill(HIST("Phi/hYResponse"), yGen, yRec);
      registry.fill(HIST("Phi/hMassVsPtGen"), ptGen, massRec);
      if (ptGen > 0.f) {
        const float ptResolution = (ptGen - ptRec) / ptGen;
        registry.fill(HIST("Phi/hPtResolution"), ptResolution);
        registry.fill(HIST("Phi/hPtGenVsResolution"), ptGen, ptResolution);
        registry.fill(HIST("Phi/hPtResponse"), ptGen, ptRec);
      }
    }
  }

  template <typename TCollision, typename TTrack>
  int getSelectedMotherId(TCollision const& collision, TTrack const& track, aod::EMMCParticles const& mcParticles,
                          std::unordered_map<int, bool> const& bestMatchMap)
  {
    if (!dimuonCut.template IsSelectedTrack<false>(track)) {
      return -1;
    }
    const auto bestMatch = bestMatchMap.find(track.globalIndex());
    if (bestMatch == bestMatchMap.end() || !bestMatch->second || track.emmcparticleId() < 0) {
      return -1;
    }

    const auto mcTrack = mcParticles.iteratorAt(track.emmcparticleId());
    if (std::abs(mcTrack.pdgCode()) != 13 || !mcTrack.has_mothers() || !isGeneratedMuonAccepted(mcTrack)) {
      return -1;
    }
    if (!(mcTrack.isPhysicalPrimary() || mcTrack.producedByGenerator())) {
      return -1;
    }
    if (cfgRequireTrueAssociation.value && mcTrack.emmceventId() != collision.emmceventId()) {
      return -1;
    }
    const auto mcCollision = mcTrack.template emmcevent_as<MyMCCollisions>();
    if (cfgEventGeneratorType.value >= 0 && mcCollision.getSubGeneratorId() != cfgEventGeneratorType.value) {
      return -1;
    }

    const int motherId = mcTrack.mothersIds()[0];
    const auto mother = mcParticles.iteratorAt(motherId);
    const int absMotherPdg = std::abs(mother.pdgCode());
    if ((absMotherPdg != 223 && absMotherPdg != 333) || !(mother.isPhysicalPrimary() || mother.producedByGenerator()) || mother.daughtersIds().size() != 2) {
      return -1;
    }
    return motherId;
  }

  template <typename TTrack1, typename TTrack2>
  void processPair(TTrack1 const& positive, TTrack2 const& negative, aod::EMMCParticles const& mcParticles)
  {
    if (!isRequestedMatchingCategory(positive, negative) || !dimuonCut.IsSelectedPair(positive, negative)) {
      return;
    }

    const auto positiveMc = mcParticles.iteratorAt(positive.emmcparticleId());
    const auto negativeMc = mcParticles.iteratorAt(negative.emmcparticleId());
    if (positiveMc.emmceventId() != negativeMc.emmceventId()) {
      return;
    }

    int motherId = aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(positiveMc, negativeMc, -13, 13, 223, mcParticles);
    if (motherId < 0) {
      motherId = aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(negativeMc, positiveMc, -13, 13, 223, mcParticles);
    }
    int motherPdg = 223;
    if (motherId < 0) {
      motherId = aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(positiveMc, negativeMc, -13, 13, 333, mcParticles);
      if (motherId < 0) {
        motherId = aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(negativeMc, positiveMc, -13, 13, 333, mcParticles);
      }
      motherPdg = 333;
    }
    if (motherId < 0) {
      return;
    }

    ROOT::Math::PtEtaPhiMVector positiveMcVector(positiveMc.pt(), positiveMc.eta(), positiveMc.phi(), constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector negativeMcVector(negativeMc.pt(), negativeMc.eta(), negativeMc.phi(), constants::physics::MassMuon);
    const auto generatedPair = positiveMcVector + negativeMcVector;
    if (generatedPair.Rapidity() < dimuoncuts.cfg_min_pair_y.value || dimuoncuts.cfg_max_pair_y.value < generatedPair.Rapidity()) {
      return;
    }

    ROOT::Math::PtEtaPhiMVector positiveVector(positive.pt(), positive.eta(), positive.phi(), constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector negativeVector(negative.pt(), negative.eta(), negative.phi(), constants::physics::MassMuon);
    const auto reconstructedPair = positiveVector + negativeVector;
    const auto mother = mcParticles.iteratorAt(motherId);
    if (motherPdg == 223) {
      fillResolutionHistograms<223>(mother.pt(), mother.y(), reconstructedPair.Pt(), reconstructedPair.Rapidity(), reconstructedPair.M());
    } else {
      fillResolutionHistograms<333>(mother.pt(), mother.y(), reconstructedPair.Pt(), reconstructedPair.Rapidity(), reconstructedPair.M());
    }
  }

  void process(FilteredMyCollisions const& collisions, MyMCCollisions const&, aod::EMMCParticles const& mcParticles, FilteredMyMCMuons const& muons)
  {
    const auto bestMatchMap = aod::pwgem::dilepton::utils::emtrackutil::findBestMatchMap(muons, dimuonCut);
    std::unordered_map<int, TracksFromEvent> tracksByEvent;
    tracksByEvent.reserve(muons.size());
    for (const auto& track : muons) {
      auto& eventTracks = tracksByEvent[track.emeventId()];
      if (track.sign() > 0) {
        eventTracks.positiveTrackIds.push_back(track.globalIndex());
      } else if (track.sign() < 0) {
        eventTracks.negativeTrackIds.push_back(track.globalIndex());
      }
    }

    for (const auto& collision : collisions) {
      const std::array<float, 3> centralities{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator.value] < eventcuts.cfgCentMin.value || eventcuts.cfgCentMax.value < centralities[cfgCentEstimator.value] || !eventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT.value && !rctChecker.checkTable(collision)) {
        continue;
      }

      const auto eventTracks = tracksByEvent.find(collision.globalIndex());
      if (eventTracks == tracksByEvent.end()) {
        continue;
      }

      std::unordered_map<int, TracksFromMother> tracksByMother;
      tracksByMother.reserve(eventTracks->second.positiveTrackIds.size() + eventTracks->second.negativeTrackIds.size());

      for (const auto trackId : eventTracks->second.positiveTrackIds) {
        const auto track = muons.rawIteratorAt(trackId);
        const int motherId = getSelectedMotherId(collision, track, mcParticles, bestMatchMap);
        if (motherId >= 0) {
          tracksByMother[motherId].positiveTrackIds.push_back(track.globalIndex());
        }
      }
      for (const auto trackId : eventTracks->second.negativeTrackIds) {
        const auto track = muons.rawIteratorAt(trackId);
        const int motherId = getSelectedMotherId(collision, track, mcParticles, bestMatchMap);
        if (motherId >= 0) {
          tracksByMother[motherId].negativeTrackIds.push_back(track.globalIndex());
        }
      }

      for (const auto& [motherId, tracks] : tracksByMother) {
        static_cast<void>(motherId); // The grouping key has already been validated above.
        for (const auto positiveId : tracks.positiveTrackIds) {
          const auto positive = muons.rawIteratorAt(positiveId);
          for (const auto negativeId : tracks.negativeTrackIds) {
            const auto negative = muons.rawIteratorAt(negativeId);
            processPair(positive, negative, mcParticles);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DimuonMcResolution>(cfgc, TaskName{"dimuon-mc-resolution"})};
}
