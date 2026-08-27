// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file globalTrackMatchingStudy.cxx
/// \brief Matching-plane residual study for reconstructed global muon tracks.

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Tasks/mftMchMatchingCategories.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TAxis.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace matchingcategory = o2::analysis::mftmchmatching;

namespace o2::aod::globaltrackmatchingstudy
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(GlobalTrackId, globalTrackId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(TrackTimeMCH, trackTimeMCH, float);
DECLARE_SOA_COLUMN(TrackTimeResMCH, trackTimeResMCH, float);
DECLARE_SOA_COLUMN(TrackTimeMFT, trackTimeMFT, float);
DECLARE_SOA_COLUMN(TrackTimeResMFT, trackTimeResMFT, float);
DECLARE_SOA_COLUMN(MatchingPlaneZ, matchingPlaneZ, float);
DECLARE_SOA_COLUMN(ZMCH, zMCH, float);
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TglMCH, tglMCH, float);
DECLARE_SOA_COLUMN(Signed1PtMCH, signed1PtMCH, float);
DECLARE_SOA_COLUMN(PtMCH, ptMCH, float);
DECLARE_SOA_COLUMN(EtaMCH, etaMCH, float);
DECLARE_SOA_COLUMN(ZMFT, zMFT, float);
DECLARE_SOA_COLUMN(XMFT, xMFT, float);
DECLARE_SOA_COLUMN(YMFT, yMFT, float);
DECLARE_SOA_COLUMN(PhiMFT, phiMFT, float);
DECLARE_SOA_COLUMN(TglMFT, tglMFT, float);
DECLARE_SOA_COLUMN(Signed1PtMFT, signed1PtMFT, float);
DECLARE_SOA_COLUMN(PtMFT, ptMFT, float);
DECLARE_SOA_COLUMN(EtaMFT, etaMFT, float);
DECLARE_SOA_COLUMN(ZGlobal, zGlobal, float);
DECLARE_SOA_COLUMN(XGlobal, xGlobal, float);
DECLARE_SOA_COLUMN(YGlobal, yGlobal, float);
DECLARE_SOA_COLUMN(PhiGlobal, phiGlobal, float);
DECLARE_SOA_COLUMN(TglGlobal, tglGlobal, float);
DECLARE_SOA_COLUMN(Signed1PtGlobal, signed1PtGlobal, float);
DECLARE_SOA_COLUMN(PtGlobal, ptGlobal, float);
DECLARE_SOA_COLUMN(EtaGlobal, etaGlobal, float);
DECLARE_SOA_COLUMN(PtTruthMCH, ptTruthMCH, float);
DECLARE_SOA_COLUMN(EtaTruthMCH, etaTruthMCH, float);
DECLARE_SOA_COLUMN(PhiTruthMCH, phiTruthMCH, float);
DECLARE_SOA_COLUMN(PtTruthMFT, ptTruthMFT, float);
DECLARE_SOA_COLUMN(EtaTruthMFT, etaTruthMFT, float);
DECLARE_SOA_COLUMN(PhiTruthMFT, phiTruthMFT, float);
DECLARE_SOA_COLUMN(PtTruthGlobal, ptTruthGlobal, float);
DECLARE_SOA_COLUMN(EtaTruthGlobal, etaTruthGlobal, float);
DECLARE_SOA_COLUMN(PhiTruthGlobal, phiTruthGlobal, float);
DECLARE_SOA_COLUMN(DeltaX, deltaX, float);
DECLARE_SOA_COLUMN(DeltaY, deltaY, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaTgl, deltaTgl, float);
DECLARE_SOA_COLUMN(DeltaSigned1Pt, deltaSigned1Pt, float);
DECLARE_SOA_COLUMN(MatchingChi2, matchingChi2, float);
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, int);
DECLARE_SOA_COLUMN(McMatchStudyCategory, mcMatchStudyCategory, int);
DECLARE_SOA_COLUMN(McPdgCodeMCH, mcPdgCodeMCH, int);
DECLARE_SOA_COLUMN(McPdgCodeMFT, mcPdgCodeMFT, int);
DECLARE_SOA_COLUMN(McFirstMotherPdgCodeMCH, mcFirstMotherPdgCodeMCH, int);
} // namespace o2::aod::globaltrackmatchingstudy

namespace o2::aod
{
DECLARE_SOA_TABLE(GlobalTrackMatchingStudies, "AOD", "GTRKMATCHSTUDY",
                  globaltrackmatchingstudy::CollisionId,
                  globaltrackmatchingstudy::GlobalTrackId,
                  globaltrackmatchingstudy::MchTrackId,
                  globaltrackmatchingstudy::MftTrackId,
                  globaltrackmatchingstudy::TrackTimeMCH,
                  globaltrackmatchingstudy::TrackTimeResMCH,
                  globaltrackmatchingstudy::TrackTimeMFT,
                  globaltrackmatchingstudy::TrackTimeResMFT,
                  globaltrackmatchingstudy::MatchingPlaneZ,
                  globaltrackmatchingstudy::ZMCH,
                  globaltrackmatchingstudy::XMCH,
                  globaltrackmatchingstudy::YMCH,
                  globaltrackmatchingstudy::PhiMCH,
                  globaltrackmatchingstudy::TglMCH,
                  globaltrackmatchingstudy::Signed1PtMCH,
                  globaltrackmatchingstudy::PtMCH,
                  globaltrackmatchingstudy::EtaMCH,
                  globaltrackmatchingstudy::ZMFT,
                  globaltrackmatchingstudy::XMFT,
                  globaltrackmatchingstudy::YMFT,
                  globaltrackmatchingstudy::PhiMFT,
                  globaltrackmatchingstudy::TglMFT,
                  globaltrackmatchingstudy::Signed1PtMFT,
                  globaltrackmatchingstudy::PtMFT,
                  globaltrackmatchingstudy::EtaMFT,
                  globaltrackmatchingstudy::ZGlobal,
                  globaltrackmatchingstudy::XGlobal,
                  globaltrackmatchingstudy::YGlobal,
                  globaltrackmatchingstudy::PhiGlobal,
                  globaltrackmatchingstudy::TglGlobal,
                  globaltrackmatchingstudy::Signed1PtGlobal,
                  globaltrackmatchingstudy::PtGlobal,
                  globaltrackmatchingstudy::EtaGlobal,
                  globaltrackmatchingstudy::PtTruthMCH,
                  globaltrackmatchingstudy::EtaTruthMCH,
                  globaltrackmatchingstudy::PhiTruthMCH,
                  globaltrackmatchingstudy::PtTruthMFT,
                  globaltrackmatchingstudy::EtaTruthMFT,
                  globaltrackmatchingstudy::PhiTruthMFT,
                  globaltrackmatchingstudy::PtTruthGlobal,
                  globaltrackmatchingstudy::EtaTruthGlobal,
                  globaltrackmatchingstudy::PhiTruthGlobal,
                  globaltrackmatchingstudy::DeltaX,
                  globaltrackmatchingstudy::DeltaY,
                  globaltrackmatchingstudy::DeltaPhi,
                  globaltrackmatchingstudy::DeltaTgl,
                  globaltrackmatchingstudy::DeltaSigned1Pt,
                  globaltrackmatchingstudy::MatchingChi2,
                  globaltrackmatchingstudy::IsCorrectMatch,
                  globaltrackmatchingstudy::McMatchStudyCategory,
                  globaltrackmatchingstudy::McPdgCodeMCH,
                  globaltrackmatchingstudy::McPdgCodeMFT,
                  globaltrackmatchingstudy::McFirstMotherPdgCodeMCH);
} // namespace o2::aod

struct GlobalTrackMatchingStudy {
  using FwdTracksData = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
  using FwdTracksMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MftTracksData = aod::MFTTracks;
  using MftTracksMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;

  Produces<aod::GlobalTrackMatchingStudies> matchingStudies;

  Configurable<float> matchingPlaneZ{"cfgMatchingPlaneZ", -77.5f, "MFT-MCH matching-plane z (cm)"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of GRPMagField"};
  Configurable<std::string> geometryPath{"geometryPath", "GLO/Config/GeometryAligned", "CCDB path of aligned geometry"};
  Configurable<bool> applyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "Read the forward-track z shift from CCDB"};
  Configurable<std::string> zShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for the forward-track z shift"};
  Configurable<float> manualZShift{"cfgManualZShift", 0.f, "Manual forward-track z shift (cm)"};
  Configurable<bool> repairMftMchRelations{"cfgRepairMftMchRelations", false,
                                           "Repair cyclic MFT/MCH component relation offsets caused by collision-unassigned suffix rows"};

  // PWGEM primary-muon skimmer selection. Dilepton analysis-level single-muon
  // cuts are intentionally not duplicated in this matching study.
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<bool> requireSel8{"cfgRequireSel8", true, "Require Sel8"};
    Configurable<float> maxAbsZvtx{"cfgMaxAbsZvtx", 10.f, "Maximum |z vertex| (cm)"};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "skimmermuoncut_group";
    Configurable<float> minPt{"minPt", 0.01f, "Minimum propagated global-muon pT (GeV/c)"};
    Configurable<float> maxPt{"maxPt", 1.e10f, "Maximum propagated global-muon pT (GeV/c)"};
    Configurable<float> minEtaGL{"minEtaGL", -1.e10f, "Minimum propagated global-muon eta"};
    Configurable<float> maxEtaGL{"maxEtaGL", 0.f, "Maximum propagated global-muon eta"};
    Configurable<float> minRabsGL{"minRabsGL", 17.6f, "Minimum R at absorber end for global muons (cm)"};
    Configurable<float> minRabs{"minRabs", 17.6f, "Minimum R at absorber end (cm)"};
    Configurable<float> midRabs{"midRabs", 26.5f, "R at absorber end separating pDCA cuts (cm)"};
    Configurable<float> maxRabs{"maxRabs", 89.5f, "Maximum R at absorber end (cm)"};
    Configurable<float> maxDCAxy{"maxDCAxy", 1.e10f, "Maximum DCAxy (cm)"};
    Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "Maximum pDCA for Rabs >= midRabs (GeV/c cm)"};
    Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "Maximum pDCA for Rabs < midRabs (GeV/c cm)"};
    Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 1.e10f, "Maximum MCH-MFT matching chi2"};
    Configurable<float> maxChi2GL{"maxChi2GL", 1.e10f, "Maximum global-muon chi2/ndf"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "Apply PWGEM MCH-momentum pT refit"};
    Configurable<float> maxDEta{"maxDEta", 1.e10f, "Maximum global-MCHMID delta eta"};
    Configurable<float> maxDPhi{"maxDPhi", 1.e10f, "Maximum global-MCHMID delta phi (rad)"};
    Configurable<int> minNmuon{"minNmuon", 0, "Minimum selected global-muon multiplicity per collision"};
  } skimmermuoncuts;

  ConfigurableAxis axisX{"axisX", {400, -20., 20.}, "#it{x} (cm)"};
  ConfigurableAxis axisY{"axisY", {400, -20., 20.}, "#it{y} (cm)"};
  ConfigurableAxis axisDeltaX{"axisDeltaX", {500, -5., 5.}, "#Delta #it{x} = #it{x}_{MFT} - #it{x}_{MCH} (cm)"};
  ConfigurableAxis axisDeltaY{"axisDeltaY", {500, -5., 5.}, "#Delta #it{y} = #it{y}_{MFT} - #it{y}_{MCH} (cm)"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {500, -0.5, 0.5}, "#Delta #it{#varphi} = #it{#varphi}_{MFT} - #it{#varphi}_{MCH} (rad)"};
  ConfigurableAxis axisDeltaTgl{"axisDeltaTgl", {500, -0.5, 0.5}, "#Delta tan(#it{#lambda})"};
  ConfigurableAxis axisDeltaSigned1Pt{"axisDeltaSigned1Pt", {500, -10., 10.}, "#Delta(#it{q}/#it{p}_{T}) (#it{c}/GeV)"};
  ConfigurableAxis axisMatchingChi2{"axisMatchingChi2", {500, 0., 100.}, "matching #chi^{2}"};
  ConfigurableAxis axisGlobalPt{"axisGlobalPt", {200, 0., 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisGlobalEta{"axisGlobalEta", {80, -4., -2.}, "#it{#eta}"};
  ConfigurableAxis axisGlobalPhi{"axisGlobalPhi", {180, 0., 2. * M_PI}, "#it{#varphi} (rad)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::globaltracking::MatchGlobalFwd matching;
  int runNumber = -1;
  float bzAtMftCenter = 0.f;
  float bzForPwgemKinematics = 0.f;
  float zShift = 0.f;

  struct TrackValues {
    float z = 0.f;
    float x = 0.f;
    float y = 0.f;
    float phi = 0.f;
    float tgl = 0.f;
    float signed1Pt = 0.f;
    float pt = 0.f;
    float eta = 0.f;
  };

  struct GlobalTrackKinematics {
    float pt = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    float dcaX = 0.f;
    float dcaY = 0.f;
    float cXXatDCA = 0.f;
    float cXYatDCA = 0.f;
    float cYYatDCA = 0.f;
    float pDCA = 0.f;
    float etaMatchedMCHMID = 0.f;
    float phiMatchedMCHMID = 0.f;
  };

  struct TruthKinematics {
    float pt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float phi = std::numeric_limits<float>::quiet_NaN();
  };

  struct McOrigin {
    bool hasParticle = false;
    int64_t particleId = -1;
    int64_t collisionId = -1;
    int pdgCode = 0;
    bool isMuon = false;
    bool isPhysicalPrimary = false;
    bool hasPiKMother = false;
    int firstMotherPdgCode = 0;
    std::unordered_set<int64_t> piKMotherIds;
  };

  struct McClassification {
    int isCorrect = -1;
    int category = matchingcategory::serialized(matchingcategory::McMatchStudyCategory::Unknown);
    int mchPdgCode = 0;
    int mftPdgCode = 0;
    int mchFirstMotherPdgCode = 0;
  };

  struct RelationLayout {
    int64_t nRows = 0;
    int64_t nUnassignedSuffixRows = 0;
    bool hasContiguousUnassignedSuffix = true;
  };

  template <typename TTracks>
  static RelationLayout inspectRelationLayout(const TTracks& tracks)
  {
    RelationLayout layout;
    layout.nRows = tracks.size();
    bool reachedUnassignedSuffix = false;
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        reachedUnassignedSuffix = true;
        ++layout.nUnassignedSuffixRows;
      } else if (reachedUnassignedSuffix) {
        layout.hasContiguousUnassignedSuffix = false;
      }
    }
    return layout;
  }

  static int64_t repairedRelationId(int64_t rawRelationId,
                                    const RelationLayout& layout,
                                    const char* relationName)
  {
    if (rawRelationId < 0 || rawRelationId >= layout.nRows) {
      LOGF(fatal, "%s relation ID %lld is outside the target table range [0, %lld)",
           relationName,
           static_cast<long long>(rawRelationId),
           static_cast<long long>(layout.nRows));
      return -1;
    }

    // The affected layout numbers relations as [unassigned][associated], while
    // the persisted target table is ordered as [associated][unassigned].
    if (rawRelationId < layout.nUnassignedSuffixRows) {
      return layout.nRows - layout.nUnassignedSuffixRows + rawRelationId;
    }
    return rawRelationId - layout.nUnassignedSuffixRows;
  }

  enum class HistogramGroup {
    All,
    CorrectPrimaryMuon,
    CorrectPiKDecayMuon,
    CorrectOtherMuon,
    FalseWithTrueCounterpart,
    NoTrueMftCounterpart
  };

  void addStudyHistograms(const std::string& prefix)
  {
    registry.add(prefix + "xMCH", "Global-track MCH component at matching plane;#it{x}_{MCH} (cm);tracks", HistType::kTH1F, {axisX});
    registry.add(prefix + "yMCH", "Global-track MCH component at matching plane;#it{y}_{MCH} (cm);tracks", HistType::kTH1F, {axisY});
    registry.add(prefix + "xMFT", "Global-track MFT component at matching plane;#it{x}_{MFT} (cm);tracks", HistType::kTH1F, {axisX});
    registry.add(prefix + "yMFT", "Global-track MFT component at matching plane;#it{y}_{MFT} (cm);tracks", HistType::kTH1F, {axisY});
    registry.add(prefix + "deltaX", "Global-track matching-plane residual;#Delta #it{x} (cm);tracks", HistType::kTH1F, {axisDeltaX});
    registry.add(prefix + "deltaY", "Global-track matching-plane residual;#Delta #it{y} (cm);tracks", HistType::kTH1F, {axisDeltaY});
    registry.add(prefix + "deltaPhi", "Global-track matching-plane residual;#Delta #it{#varphi} (rad);tracks", HistType::kTH1F, {axisDeltaPhi});
    registry.add(prefix + "deltaTgl", "Global-track matching-plane residual;#Delta tan(#it{#lambda});tracks", HistType::kTH1F, {axisDeltaTgl});
    registry.add(prefix + "deltaSigned1Pt", "Global-track matching-plane residual;#Delta(#it{q}/#it{p}_{T}) (#it{c}/GeV);tracks", HistType::kTH1F, {axisDeltaSigned1Pt});
    registry.add(prefix + "matchingChi2", "Global-track MFT-MCH matching #chi^{2};matching #chi^{2};tracks", HistType::kTH1F, {axisMatchingChi2});
    registry.add(prefix + "deltaXVsMatchingChi2", "Matching-plane residual versus matching #chi^{2};matching #chi^{2};#Delta #it{x} (cm)", HistType::kTH2F, {axisMatchingChi2, axisDeltaX});
    registry.add(prefix + "deltaYVsMatchingChi2", "Matching-plane residual versus matching #chi^{2};matching #chi^{2};#Delta #it{y} (cm)", HistType::kTH2F, {axisMatchingChi2, axisDeltaY});
    registry.add(prefix + "globalPt", "Global-track kinematics at the primary vertex;#it{p}_{T} (GeV/#it{c});tracks", HistType::kTH1F, {axisGlobalPt});
    registry.add(prefix + "globalEta", "Global-track kinematics at the primary vertex;#it{#eta};tracks", HistType::kTH1F, {axisGlobalEta});
    registry.add(prefix + "globalPhi", "Global-track kinematics at the primary vertex;#it{#varphi} (rad);tracks", HistType::kTH1F, {axisGlobalPhi});
    registry.add(prefix + "globalEtaPhi", "Global-track kinematics at the primary vertex;#it{#varphi} (rad);#it{#eta}", HistType::kTH2F, {axisGlobalPhi, axisGlobalEta});
  }

  void init(InitContext const&)
  {
    const int nEnabledProcesses = static_cast<int>(doprocessData) + static_cast<int>(doprocessDataWithFTTCA) +
                                  static_cast<int>(doprocessMC) + static_cast<int>(doprocessMCWithFTTCA);
    if (nEnabledProcesses != 1) {
      LOGF(fatal, "Enable exactly one data/MC process, with or without FwdTrackAssoc");
    }
    if (!std::isfinite(matchingPlaneZ.value)) {
      LOGF(fatal, "cfgMatchingPlaneZ must be finite");
    }
    if (!std::isfinite(manualZShift.value) || !std::isfinite(eventcuts.maxAbsZvtx.value) || eventcuts.maxAbsZvtx.value < 0.f ||
        !std::isfinite(skimmermuoncuts.minPt.value) || !std::isfinite(skimmermuoncuts.maxPt.value) ||
        skimmermuoncuts.minPt.value > skimmermuoncuts.maxPt.value ||
        !std::isfinite(skimmermuoncuts.minEtaGL.value) || !std::isfinite(skimmermuoncuts.maxEtaGL.value) ||
        skimmermuoncuts.minEtaGL.value > skimmermuoncuts.maxEtaGL.value ||
        !std::isfinite(skimmermuoncuts.minRabs.value) || !std::isfinite(skimmermuoncuts.midRabs.value) ||
        !std::isfinite(skimmermuoncuts.maxRabs.value) || skimmermuoncuts.minRabs.value > skimmermuoncuts.midRabs.value ||
        skimmermuoncuts.midRabs.value > skimmermuoncuts.maxRabs.value ||
        !std::isfinite(skimmermuoncuts.maxDCAxy.value) || skimmermuoncuts.maxDCAxy.value < 0.f ||
        !std::isfinite(skimmermuoncuts.maxDEta.value) || skimmermuoncuts.maxDEta.value <= 0.f ||
        !std::isfinite(skimmermuoncuts.maxDPhi.value) || skimmermuoncuts.maxDPhi.value <= 0.f ||
        skimmermuoncuts.minNmuon.value < 0) {
      LOGF(fatal, "Invalid PWGEM-style event or skimmer-muon cut configuration");
    }
    if (repairMftMchRelations.value) {
      LOGF(warning,
           "cfgRepairMftMchRelations is enabled: applying the opt-in cyclic MFT/MCH relation repair "
           "after validating the collision-unassigned suffix layout");
    }

    ccdbManager->setURL(ccdbUrl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();

    addStudyHistograms("all/");
    if (doprocessDataWithFTTCA || doprocessMCWithFTTCA) {
      auto associationStatus = registry.add<TH1>("association/type0Status",
                                                  "Type-0 GlobalMuonTrack collision-association status;status;tracks",
                                                  HistType::kTH1I, {{2, -0.5, 1.5}});
      associationStatus->GetXaxis()->SetBinLabel(1, "all type-0 tracks");
      associationStatus->GetXaxis()->SetBinLabel(2, ">= 1 FwdTrackAssoc");
      registry.add("association/nAssociationsPerType0Track",
                   "FwdTrackAssoc multiplicity per type-0 GlobalMuonTrack;collision associations;tracks",
                   HistType::kTH1I, {{21, -0.5, 20.5}});
    }
    if (doprocessMC || doprocessMCWithFTTCA) {
      addStudyHistograms("mc/correctPrimaryMuon/");
      addStudyHistograms("mc/correctPiKDecayMuon/");
      addStudyHistograms("mc/correctOtherMuon/");
      addStudyHistograms("mc/falseWithTrueCounterpart/");
      addStudyHistograms("mc/noTrueMftCounterpart/");
      auto categoryHistogram = registry.add<TH1>("mc/category", "Global-track MC matching category;category;tracks", HistType::kTH1F, {{5, -0.5, 4.5}});
      registry.add("mc/nOutputCandidatesPerMchMuon",
                   "Output GlobalMuonTrack candidates per MC MCH-side muon;number of output candidates per MC muon;MC muons",
                   HistType::kTH1I, {{51, -0.5, 50.5}});
      categoryHistogram->GetXaxis()->SetBinLabel(1, "correct primary muon");
      categoryHistogram->GetXaxis()->SetBinLabel(2, "correct #pi/K decay muon");
      categoryHistogram->GetXaxis()->SetBinLabel(3, "correct other muon");
      categoryHistogram->GetXaxis()->SetBinLabel(4, "false with true counterpart");
      categoryHistogram->GetXaxis()->SetBinLabel(5, "truth muon without reconstructed MFT counterpart");
    }
  }

  template <typename TBCs>
  void prepareRun(const TBCs& bcs)
  {
    if (bcs.size() == 0) {
      LOGF(fatal, "Cannot propagate global-track components without BC timestamps for the magnetic field");
    }
    const auto bc = bcs.begin();
    if (runNumber == bc.runNumber()) {
      return;
    }
    auto* grpmag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    if (grpmag == nullptr) {
      LOGF(fatal, "Could not load GRPMagField from %s for timestamp %lld", grpmagPath.value.c_str(), static_cast<long long>(bc.timestamp()));
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Loading aligned geometry from CCDB path %s", geometryPath.value.c_str());
      auto* geometry = ccdbManager->getForTimeStamp<TGeoManager>(geometryPath, bc.timestamp());
      if (geometry == nullptr || !o2::base::GeometryManager::isGeometryLoaded()) {
        LOGF(fatal, "Could not load aligned geometry from %s for timestamp %lld", geometryPath.value.c_str(), static_cast<long long>(bc.timestamp()));
      }
    }
    o2::mch::TrackExtrap::setField();
    bzAtMftCenter = grpmag->getNominalL3Field();
    constexpr double centerMFT[3] = {0., 0., -61.4};
    auto* magneticField = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (magneticField == nullptr) {
      LOGF(fatal, "Magnetic field is not initialized for PWGEM-style global-track propagation");
    }
    bzForPwgemKinematics = magneticField->getBz(centerMFT);
    if (applyZShiftFromCCDB) {
      const auto* ccdbZShift = ccdbManager->getForTimeStamp<std::vector<float>>(zShiftPath, bc.timestamp());
      zShift = ccdbZShift != nullptr && !ccdbZShift->empty() ? ccdbZShift->front() : 0.f;
    } else {
      zShift = manualZShift;
    }
    runNumber = bc.runNumber();
  }

  static TrackValues values(const o2::track::TrackParCovFwd& track)
  {
    TrackValues result;
    result.z = track.getZ();
    result.x = track.getX();
    result.y = track.getY();
    result.phi = track.getPhi();
    result.tgl = track.getTgl();
    result.signed1Pt = track.getInvQPt();
    constexpr float minimumAbsSigned1Pt = 1.e-12f;
    result.pt = std::abs(result.signed1Pt) > minimumAbsSigned1Pt ? std::abs(1.f / result.signed1Pt) : 0.f;
    result.eta = std::asinh(result.tgl);
    return result;
  }

  template <typename TTrack>
  static TruthKinematics truthKinematics(const TTrack& track)
  {
    if (!track.has_mcParticle()) {
      return {};
    }
    const auto& particle = track.mcParticle();
    return {particle.pt(), particle.eta(), particle.phi()};
  }

  template <typename TGlobalTrack, typename TMchTrack, typename TCollision>
  std::optional<GlobalTrackKinematics> pwgemGlobalTrackKinematics(const TGlobalTrack& globalTrack,
                                                                    const TMchTrack& mchTrack,
                                                                    const TCollision& collision) const
  {
    // Match the active PWGEM primary-muon skimmer convention: propagate the
    // global track to the PV, then use the matched MCH momentum for the default
    // global-muon pT refit while retaining the propagated global eta and phi.
    const auto globalAtPV = o2::aod::fwdtrackutils::propagateMuon(
      globalTrack, globalTrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex,
      matchingPlaneZ.value, bzForPwgemKinematics, zShift);
    const auto mchAtPV = o2::aod::fwdtrackutils::propagateMuon(
      mchTrack, mchTrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex,
      matchingPlaneZ.value, bzForPwgemKinematics, zShift);

    GlobalTrackKinematics result;
    result.eta = globalAtPV.getEta();
    result.phi = globalAtPV.getPhi();
    o2::math_utils::bringTo02Pi(result.phi);
    result.etaMatchedMCHMID = mchAtPV.getEta();
    result.phiMatchedMCHMID = mchAtPV.getPhi();
    o2::math_utils::bringTo02Pi(result.phiMatchedMCHMID);
    result.pt = globalAtPV.getPt();

    auto globalAtDCA = o2::aod::fwdtrackutils::propagateMuon(
      globalTrack, globalTrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA,
      matchingPlaneZ.value, bzForPwgemKinematics, zShift);
    if (skimmermuoncuts.refitGlobalMuon.value) {
      if (!std::isfinite(globalAtPV.getPt()) || globalAtPV.getPt() <= 0.f) {
        return std::nullopt;
      }
      const float refittedPt = mchAtPV.getP() * std::sin(2.f * std::atan(std::exp(-result.eta)));
      if (!std::isfinite(refittedPt) || refittedPt <= 0.f) {
        return std::nullopt;
      }
      const float sfPt = refittedPt / globalAtPV.getPt();
      const auto& covarianceAtPV = globalAtPV.getCovariances();
      const auto refittedGlobalAtPV = o2::aod::fwdtrackutils::getTrackParCovFwdShiftManual(
        globalAtPV.getX(), globalAtPV.getY(), globalAtPV.getPhi(), globalAtPV.getTgl(), globalAtPV.getInvQPt() / sfPt,
        covarianceAtPV(0, 0),
        covarianceAtPV(1, 0), covarianceAtPV(1, 1),
        covarianceAtPV(2, 0), covarianceAtPV(2, 1), covarianceAtPV(2, 2),
        covarianceAtPV(3, 0), covarianceAtPV(3, 1), covarianceAtPV(3, 2), covarianceAtPV(3, 3),
        covarianceAtPV(4, 0) / sfPt, covarianceAtPV(4, 1) / sfPt, covarianceAtPV(4, 2) / sfPt,
        covarianceAtPV(4, 3) / sfPt, covarianceAtPV(4, 4) / sfPt / sfPt,
        globalAtPV.getZ(), 0.f, globalTrack.chi2());
      globalAtDCA = o2::aod::fwdtrackutils::propagateTrackParCovFwd(
        refittedGlobalAtPV, globalTrack.trackType(), collision,
        o2::aod::fwdtrackutils::propagationPoint::kToDCA, matchingPlaneZ.value, bzForPwgemKinematics);
      result.pt = refittedPt;
    }

    const auto mchAtDCA = o2::aod::fwdtrackutils::propagateMuon(
      mchTrack, mchTrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA,
      matchingPlaneZ.value, bzForPwgemKinematics, zShift);
    const float mchDcaX = mchAtDCA.getX() - collision.posX();
    const float mchDcaY = mchAtDCA.getY() - collision.posY();
    result.pDCA = mchTrack.p() * std::sqrt(mchDcaX * mchDcaX + mchDcaY * mchDcaY);
    result.dcaX = globalAtDCA.getX() - collision.posX();
    result.dcaY = globalAtDCA.getY() - collision.posY();
    const auto& covarianceAtDCA = globalAtDCA.getCovariances();
    result.cXXatDCA = covarianceAtDCA(0, 0);
    result.cXYatDCA = covarianceAtDCA(1, 0);
    result.cYYatDCA = covarianceAtDCA(1, 1);

    if (!std::isfinite(globalAtDCA.getPhi()) || !std::isfinite(globalAtDCA.getTgl()) || std::abs(globalAtDCA.getTgl()) < 1.e-8f) {
      return std::nullopt;
    }
    const float tx = std::cos(globalAtDCA.getPhi()) / globalAtDCA.getTgl();
    const float ty = std::sin(globalAtDCA.getPhi()) / globalAtDCA.getTgl();
    result.cXXatDCA += collision.covXX() - 2.f * tx * collision.covXZ() + tx * tx * collision.covZZ();
    result.cXYatDCA += collision.covXY() - ty * collision.covXZ() - tx * collision.covYZ() + tx * ty * collision.covZZ();
    result.cYYatDCA += collision.covYY() - 2.f * ty * collision.covYZ() + ty * ty * collision.covZZ();
    if (!std::isfinite(result.dcaX) || !std::isfinite(result.dcaY) || !std::isfinite(result.pDCA) ||
        !std::isfinite(result.cXXatDCA) || !std::isfinite(result.cXYatDCA) || !std::isfinite(result.cYYatDCA)) {
      return std::nullopt;
    }
    return result;
  }

  template <typename TCollision>
  bool isSelectedEvent(const TCollision& collision) const
  {
    return (!eventcuts.requireSel8.value || collision.sel8()) && std::abs(collision.posZ()) < eventcuts.maxAbsZvtx.value;
  }

  template <typename TGlobalTrack, typename TMftTrack>
  bool passesSkimmerMuonSelection(const TGlobalTrack& globalTrack, const TMftTrack& mftTrack, const GlobalTrackKinematics& global) const
  {
    const float rabs = globalTrack.rAtAbsorberEnd();
    const float dcaXY = std::sqrt(global.dcaX * global.dcaX + global.dcaY * global.dcaY);
    const float ndf = 2.f * (globalTrack.nClusters() + mftTrack.nClusters()) - 5.f;
    const float globalChi2Ndf = globalTrack.chi2() / ndf;
    if (ndf <= 0.f || !std::isfinite(rabs) || !std::isfinite(dcaXY) || !std::isfinite(globalChi2Ndf) ||
        globalTrack.chi2MatchMCHMID() < 0.f || globalTrack.chi2MatchMCHMFT() < 0.f || globalTrack.chi2() < 0.f) {
      return false;
    }
    if (global.pt < skimmermuoncuts.minPt.value || skimmermuoncuts.maxPt.value < global.pt ||
        global.eta < skimmermuoncuts.minEtaGL.value || skimmermuoncuts.maxEtaGL.value < global.eta ||
        rabs < skimmermuoncuts.minRabs.value || skimmermuoncuts.maxRabs.value < rabs ||
        rabs < skimmermuoncuts.minRabsGL.value || skimmermuoncuts.maxRabs.value < rabs ||
        global.pDCA > (rabs < skimmermuoncuts.midRabs.value ? skimmermuoncuts.maxPDCAforSmallR.value : skimmermuoncuts.maxPDCAforLargeR.value) ||
        dcaXY > skimmermuoncuts.maxDCAxy.value || globalChi2Ndf > skimmermuoncuts.maxChi2GL.value ||
        globalTrack.chi2MatchMCHMFT() > skimmermuoncuts.maxMatchingChi2MCHMFT.value) {
      return false;
    }
    float dphi = global.phiMatchedMCHMID - global.phi;
    o2::math_utils::bringToPMPi(dphi);
    const float deta = global.etaMatchedMCHMID - global.eta;
    if (std::sqrt(std::pow(deta / skimmermuoncuts.maxDEta.value, 2) + std::pow(dphi / skimmermuoncuts.maxDPhi.value, 2)) > 1.f) {
      return false;
    }
    return true;
  }

  std::optional<o2::track::TrackParCovFwd> propagateMchToMatchingPlane(const o2::track::TrackParCovFwd& track)
  {
    o2::dataformats::GlobalFwdTrack globalTrack;
    globalTrack.setParameters(track.getParameters());
    globalTrack.setZ(track.getZ());
    globalTrack.setCovariances(track.getCovariances());
    auto mchTrack = matching.FwdtoMCH(globalTrack);

    constexpr float absorberFrontZ = -90.f;
    const bool propagated = matchingPlaneZ.value >= absorberFrontZ ? o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, matchingPlaneZ.value)
                                                                   : o2::mch::TrackExtrap::extrapToZCov(mchTrack, matchingPlaneZ.value);
    if (!propagated) {
      return std::nullopt;
    }
    return matching.MCHtoFwd(mchTrack);
  }

  o2::track::TrackParCovFwd propagateMftToMatchingPlane(const o2::track::TrackParCovFwd& track) const
  {
    auto propagated = track;
    propagated.propagateToZ(matchingPlaneZ.value, bzAtMftCenter);
    return propagated;
  }

  static float wrappedDeltaPhi(float mftPhi, float mchPhi)
  {
    return std::atan2(std::sin(mftPhi - mchPhi), std::cos(mftPhi - mchPhi));
  }

  template <typename TTrack>
  static McOrigin mcOrigin(const TTrack& track)
  {
    McOrigin origin;
    if (!track.has_mcParticle()) {
      return origin;
    }
    const auto& particle = track.mcParticle();
    origin.hasParticle = true;
    origin.particleId = particle.globalIndex();
    origin.collisionId = particle.mcCollisionId();
    origin.pdgCode = particle.pdgCode();
    origin.isMuon = std::abs(particle.pdgCode()) == 13;
    origin.isPhysicalPrimary = particle.isPhysicalPrimary();
    const auto& mothers = particle.template mothers_as<aod::McParticles>();
    if (!mothers.empty()) {
      origin.firstMotherPdgCode = mothers[0].pdgCode();
    }
    for (const auto& mother : mothers) {
      const int absPdgCode = std::abs(mother.pdgCode());
      if (absPdgCode == 211 || absPdgCode == 321) {
        origin.hasPiKMother = true;
        origin.piKMotherIds.insert(mother.globalIndex());
      }
    }
    return origin;
  }

  static bool containsParticle(const std::unordered_map<int64_t, std::unordered_set<int64_t>>& particlesByCollision,
                               int64_t collisionId,
                               int64_t particleId)
  {
    const auto collision = particlesByCollision.find(collisionId);
    return collisionId >= 0 && particleId >= 0 && collision != particlesByCollision.end() && collision->second.contains(particleId);
  }

  static bool containsAnyParticle(const std::unordered_map<int64_t, std::unordered_set<int64_t>>& particlesByCollision,
                                  int64_t collisionId,
                                  const std::unordered_set<int64_t>& particleIds)
  {
    const auto collision = particlesByCollision.find(collisionId);
    if (collisionId < 0 || collision == particlesByCollision.end()) {
      return false;
    }
    for (const auto particleId : particleIds) {
      if (collision->second.contains(particleId)) {
        return true;
      }
    }
    return false;
  }

  static McClassification classify(const McOrigin& mchOrigin,
                                   const McOrigin& mftOrigin,
                                   const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftParticlesByMcCollision)
  {
    McClassification result;
    result.mchPdgCode = mchOrigin.pdgCode;
    result.mftPdgCode = mftOrigin.pdgCode;
    result.mchFirstMotherPdgCode = mchOrigin.firstMotherPdgCode;

    // Matching-study categories are defined only for a truth-labelled MCH-side muon.
    // Keep non-muons and tracks without an MCH truth label as Unknown.
    if (!mchOrigin.hasParticle || !mchOrigin.isMuon) {
      return result;
    }

    const bool sameMcCollision = mchOrigin.hasParticle && mftOrigin.hasParticle &&
                                 mchOrigin.collisionId >= 0 && mchOrigin.collisionId == mftOrigin.collisionId;
    const bool sameMcParticle = mchOrigin.hasParticle && mftOrigin.hasParticle && mchOrigin.particleId == mftOrigin.particleId;
    const bool piKMotherMatch = mchOrigin.isMuon && mchOrigin.hasPiKMother && sameMcCollision && mchOrigin.piKMotherIds.contains(mftOrigin.particleId);
    const bool correct = matchingcategory::isCorrectPair(sameMcParticle, piKMotherMatch);
    const bool hasTrueMftCounterpart = containsParticle(mftParticlesByMcCollision, mchOrigin.collisionId, mchOrigin.particleId);
    const bool hasMotherMftCounterpart = mchOrigin.isMuon && mchOrigin.hasPiKMother &&
                                         containsAnyParticle(mftParticlesByMcCollision, mchOrigin.collisionId, mchOrigin.piKMotherIds);

    result.isCorrect = correct ? 1 : 0;
    result.category = matchingcategory::serialized(
      matchingcategory::classifyStudyPair(correct,
                                          mchOrigin.isMuon,
                                          mchOrigin.hasPiKMother,
                                          mchOrigin.isPhysicalPrimary,
                                          hasTrueMftCounterpart || hasMotherMftCounterpart));
    return result;
  }

  void fillHistograms(HistogramGroup group,
                      const TrackValues& mch,
                      const TrackValues& mft,
                      const GlobalTrackKinematics& global,
                      float deltaX,
                      float deltaY,
                      float deltaPhi,
                      float deltaTgl,
                      float deltaSigned1Pt,
                      float matchingChi2)
  {
#define FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS(prefix)                                         \
  registry.fill(HIST(prefix "xMCH"), mch.x);                                                 \
  registry.fill(HIST(prefix "yMCH"), mch.y);                                                 \
  registry.fill(HIST(prefix "xMFT"), mft.x);                                                 \
  registry.fill(HIST(prefix "yMFT"), mft.y);                                                 \
  registry.fill(HIST(prefix "deltaX"), deltaX);                                              \
  registry.fill(HIST(prefix "deltaY"), deltaY);                                              \
  registry.fill(HIST(prefix "deltaPhi"), deltaPhi);                                          \
  registry.fill(HIST(prefix "deltaTgl"), deltaTgl);                                          \
  registry.fill(HIST(prefix "deltaSigned1Pt"), deltaSigned1Pt);                              \
  registry.fill(HIST(prefix "matchingChi2"), matchingChi2);                                  \
  registry.fill(HIST(prefix "deltaXVsMatchingChi2"), matchingChi2, deltaX);                   \
  registry.fill(HIST(prefix "deltaYVsMatchingChi2"), matchingChi2, deltaY);                   \
  registry.fill(HIST(prefix "globalPt"), global.pt);                                         \
  registry.fill(HIST(prefix "globalEta"), global.eta);                                       \
  registry.fill(HIST(prefix "globalPhi"), global.phi);                                       \
  registry.fill(HIST(prefix "globalEtaPhi"), global.phi, global.eta)

    switch (group) {
      case HistogramGroup::All:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("all/");
        break;
      case HistogramGroup::CorrectPrimaryMuon:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("mc/correctPrimaryMuon/");
        break;
      case HistogramGroup::CorrectPiKDecayMuon:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("mc/correctPiKDecayMuon/");
        break;
      case HistogramGroup::CorrectOtherMuon:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("mc/correctOtherMuon/");
        break;
      case HistogramGroup::FalseWithTrueCounterpart:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("mc/falseWithTrueCounterpart/");
        break;
      case HistogramGroup::NoTrueMftCounterpart:
        FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS("mc/noTrueMftCounterpart/");
        break;
    }

#undef FILL_GLOBAL_TRACK_MATCHING_HISTOGRAMS
  }

  void fillMcCategoryHistograms(int category,
                                const TrackValues& mch,
                                const TrackValues& mft,
                                const GlobalTrackKinematics& global,
                                float deltaX,
                                float deltaY,
                                float deltaPhi,
                                float deltaTgl,
                                float deltaSigned1Pt,
                                float matchingChi2)
  {
    using Category = matchingcategory::McMatchStudyCategory;
    registry.fill(HIST("mc/category"), category);
    switch (static_cast<Category>(category)) {
      case Category::CorrectPrimaryMuon:
        fillHistograms(HistogramGroup::CorrectPrimaryMuon, mch, mft, global, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::CorrectPiKDecayMuon:
        fillHistograms(HistogramGroup::CorrectPiKDecayMuon, mch, mft, global, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::CorrectOtherMuon:
        fillHistograms(HistogramGroup::CorrectOtherMuon, mch, mft, global, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::FalseWithTrueCounterpart:
        fillHistograms(HistogramGroup::FalseWithTrueCounterpart, mch, mft, global, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::NoTrueMftCounterpart:
        fillHistograms(HistogramGroup::NoTrueMftCounterpart, mch, mft, global, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      default:
        break;
    }
  }

  using FwdTrackToCollisions = std::unordered_multimap<int64_t, int64_t>;

  static FwdTrackToCollisions makeFwdTrackToCollisions(const aod::FwdTrackAssoc& associations)
  {
    FwdTrackToCollisions result;
    result.reserve(associations.size());
    for (const auto& association : associations) {
      result.emplace(association.fwdtrackId(), association.collisionId());
    }
    return result;
  }

  template <bool isMC, typename TCollisions, typename TFwdTracks, typename TMftTracks>
  void runStudy(const TFwdTracks& fwdTracks,
                const TMftTracks& mftTracks,
                const TCollisions& collisions,
                const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftParticlesByMcCollision,
                const FwdTrackToCollisions* fwdTrackToCollisions = nullptr)
  {
    struct BestGlobalMuonCandidate {
      int64_t globalTrackId;
      float matchingChi2;
    };
    std::unordered_map<int64_t, BestGlobalMuonCandidate> bestGlobalCandidateByMchTrack;
    std::unordered_map<int64_t, int> outputCandidatesPerMcMuon;

    // A contiguous unassigned suffix is required for this repair, but does not
    // by itself prove that a production has shifted relations. Keep it opt-in.
    const bool repairRelations = repairMftMchRelations.value;
    RelationLayout fwdRelationLayout;
    RelationLayout mftRelationLayout;
    if (repairRelations) {
      fwdRelationLayout = inspectRelationLayout(fwdTracks);
      mftRelationLayout = inspectRelationLayout(mftTracks);
      if (!fwdRelationLayout.hasContiguousUnassignedSuffix || !mftRelationLayout.hasContiguousUnassignedSuffix) {
        LOGF(fatal,
             "Cannot repair MFT/MCH relations: collision-unassigned rows are not a contiguous suffix "
             "(FwdTracks=%s, MFTTracks=%s)",
             fwdRelationLayout.hasContiguousUnassignedSuffix ? "valid" : "invalid",
             mftRelationLayout.hasContiguousUnassignedSuffix ? "valid" : "invalid");
      }
    }

    auto resolveMchTrack = [&](const auto& globalTrack) {
      if (!repairRelations) {
        return globalTrack.template matchMCHTrack_as<TFwdTracks>();
      }
      const int64_t correctedId = repairedRelationId(globalTrack.matchMCHTrackId(), fwdRelationLayout, "MCH");
      return fwdTracks.rawIteratorAt(correctedId);
    };
    auto resolveMftTrack = [&](const auto& globalTrack) {
      if (!repairRelations) {
        return globalTrack.template matchMFTTrack_as<TMftTracks>();
      }
      const int64_t correctedId = repairedRelationId(globalTrack.matchMFTTrackId(), mftRelationLayout, "MFT");
      return mftTracks.rawIteratorAt(correctedId);
    };

    auto selectedCollisionIds = [&](const auto& globalTrack) {
      std::vector<int64_t> result;
      if (fwdTrackToCollisions != nullptr) {
        const auto [begin, end] = fwdTrackToCollisions->equal_range(globalTrack.globalIndex());
        for (auto association = begin; association != end; ++association) {
          const auto collision = collisions.rawIteratorAt(association->second);
          if (isSelectedEvent(collision)) {
            result.emplace_back(association->second);
          }
        }
      } else if (globalTrack.has_collision()) {
        const auto collision = collisions.rawIteratorAt(globalTrack.collisionId());
        if (isSelectedEvent(collision)) {
          result.emplace_back(globalTrack.collisionId());
        }
      }
      return result;
    };

    // Keep one GlobalMuonTrack candidate per MCH-MID relation, as done in the PWGDQ skimmers.
    // Only type 3 is a MuonStandaloneTrack (MCH-MID) suitable for MCH propagation.
    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        continue;
      }
      if (fwdTrackToCollisions != nullptr) {
        registry.fill(HIST("association/type0Status"), 0.);
        const int nAssociations = fwdTrackToCollisions->count(globalTrack.globalIndex());
        registry.fill(HIST("association/nAssociationsPerType0Track"), nAssociations);
        if (nAssociations == 0) {
          continue;
        }
        registry.fill(HIST("association/type0Status"), 1.);
      }

      if (selectedCollisionIds(globalTrack).empty()) {
        continue;
      }

      const auto mchTrack = resolveMchTrack(globalTrack);
      const int mchTrackType = static_cast<int>(mchTrack.trackType());
      if (mchTrackType != aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        // LOGF(info, "%s Skipping global track %lld: MCH relation does not resolve to type-3 MuonStandaloneTrack",
        //      isMC ? "MC" : "Real-data", static_cast<long long>(globalTrack.globalIndex()));
        continue;
      }
      const float matchingChi2 = globalTrack.chi2MatchMCHMFT();
      if (!std::isfinite(matchingChi2) || matchingChi2 <= 0.f) {
        // LOGF(info, "%s Skipping global track %lld: invalid chi2MatchMCHMFT=%.6g",
        //      isMC ? "MC" : "Real-data", static_cast<long long>(globalTrack.globalIndex()), static_cast<double>(matchingChi2));
        continue;
      }

      const int64_t mchTrackId = mchTrack.globalIndex();
      const auto bestCandidate = bestGlobalCandidateByMchTrack.find(mchTrackId);
      if (bestCandidate == bestGlobalCandidateByMchTrack.end() || matchingChi2 < bestCandidate->second.matchingChi2) {
        bestGlobalCandidateByMchTrack[mchTrackId] = {globalTrack.globalIndex(), matchingChi2};
      }
    }

    // Mirror the skimmer's per-collision multiplicity requirement before
    // materialising any selected global-muon candidates.
    std::unordered_map<int64_t, int> selectedGlobalMuonsPerCollision;
    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        continue;
      }
      const auto mchTrack = resolveMchTrack(globalTrack);
      const auto mftTrack = resolveMftTrack(globalTrack);
      const auto bestCandidate = bestGlobalCandidateByMchTrack.find(mchTrack.globalIndex());
      if (bestCandidate == bestGlobalCandidateByMchTrack.end() || bestCandidate->second.globalTrackId != globalTrack.globalIndex()) {
        continue;
      }
      for (const auto collisionId : selectedCollisionIds(globalTrack)) {
        const auto collision = collisions.rawIteratorAt(collisionId);
        const auto globalKinematics = pwgemGlobalTrackKinematics(globalTrack, mchTrack, collision);
        if (globalKinematics && passesSkimmerMuonSelection(globalTrack, mftTrack, *globalKinematics)) {
          ++selectedGlobalMuonsPerCollision[collisionId];
        }
      }
    }

    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        continue;
      }
      const auto collisionIds = selectedCollisionIds(globalTrack);
      if (collisionIds.empty()) {
        continue;
      }

      const auto mchTrack = resolveMchTrack(globalTrack);
      const auto mftTrack = resolveMftTrack(globalTrack);
      const int64_t mchTrackId = mchTrack.globalIndex();
      const int64_t mftTrackId = mftTrack.globalIndex();
      const int mchTrackType = static_cast<int>(mchTrack.trackType());
      if (mchTrackType != aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        continue;
      }
      const float matchingChi2 = globalTrack.chi2MatchMCHMFT();
      if (!std::isfinite(matchingChi2) || matchingChi2 <= 0.f) {
        continue;
      }
      const auto bestCandidate = bestGlobalCandidateByMchTrack.find(mchTrackId);
      if (bestCandidate == bestGlobalCandidateByMchTrack.end() || bestCandidate->second.globalTrackId != globalTrack.globalIndex()) {
        // LOGF(info, "%s Skipping non-best GlobalMuon candidate: global=%lld, MCH=%lld, chi2MatchMCHMFT=%.6g, best global=%lld, best chi2=%.6g",
        //      isMC ? "MC" : "Real-data", static_cast<long long>(globalTrack.globalIndex()),
        //      static_cast<long long>(mchTrackId), static_cast<double>(matchingChi2),
        //      static_cast<long long>(bestCandidate->second.globalTrackId), static_cast<double>(bestCandidate->second.matchingChi2));
        continue;
      }

      // MFT covariance rows are sparse and indexed by matchMFTTrackId. The coordinates
      // and helix propagation used here do not require them, so avoid an invalid row-wise join.
      const auto mftRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(mftTrack, 0.f);
      const auto mchRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(mchTrack, 0.f, mchTrack);
      const auto globalRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(globalTrack, 0.f, globalTrack);
      // LOGF(info, "%s MCH input: global=%lld, MCH=%lld, MFT=%lld, MCH type=%d, z=%.4f cm",
      //      isMC ? "MC" : "Real-data", static_cast<long long>(globalTrack.globalIndex()),
      //      static_cast<long long>(mchTrackId), static_cast<long long>(mftTrackId), mchTrackType,
      //      static_cast<double>(mchRaw.getZ()));
      const auto mchAtMatching = propagateMchToMatchingPlane(mchRaw);
      if (!mchAtMatching) {
        LOGF(warning, "MCH propagation failed for global track %lld", static_cast<long long>(globalTrack.globalIndex()));
        continue;
      }
      const auto mftAtMatching = propagateMftToMatchingPlane(mftRaw);
      const auto globalAtMatching = propagateMftToMatchingPlane(globalRaw);
      const auto mch = values(*mchAtMatching);
      const auto mft = values(mftAtMatching);
      const auto globalAtMatchingValues = values(globalAtMatching);

      const float deltaX = mft.x - mch.x;
      const float deltaY = mft.y - mch.y;
      const float deltaPhi = wrappedDeltaPhi(mft.phi, mch.phi);
      const float deltaTgl = mft.tgl - mch.tgl;
      const float deltaSigned1Pt = mft.signed1Pt - mch.signed1Pt;
      // This is the matching chi2 stored during reconstruction. It is not recalculated
      // when cfgMatchingPlaneZ differs from the plane used by the reconstruction.
      for (const auto collisionId : collisionIds) {
        if (selectedGlobalMuonsPerCollision[collisionId] < skimmermuoncuts.minNmuon.value) {
          continue;
        }
        const auto collision = collisions.rawIteratorAt(collisionId);
        const auto globalKinematics = pwgemGlobalTrackKinematics(globalTrack, mchTrack, collision);
        if (!globalKinematics || !passesSkimmerMuonSelection(globalTrack, mftTrack, *globalKinematics)) {
          continue;
        }
        McClassification mc;
        TruthKinematics mchTruth;
        TruthKinematics mftTruth;
        TruthKinematics globalTruth;
        if constexpr (isMC) {
          const auto mchOrigin = mcOrigin(mchTrack);
          const auto mftOrigin = mcOrigin(mftTrack);
          mchTruth = truthKinematics(mchTrack);
          mftTruth = truthKinematics(mftTrack);
          globalTruth = truthKinematics(globalTrack);
          mc = classify(mchOrigin, mftOrigin, mftParticlesByMcCollision);
          if (mchOrigin.hasParticle && mchOrigin.isMuon) {
            ++outputCandidatesPerMcMuon[mchOrigin.particleId];
          }
        }

        matchingStudies(collisionId,
                      globalTrack.globalIndex(),
                      mchTrackId,
                      mftTrackId,
                      mchTrack.trackTime(),
                      mchTrack.trackTimeRes(),
                      mftTrack.trackTime(),
                      mftTrack.trackTimeRes(),
                      matchingPlaneZ.value,
                      mch.z,
                      mch.x,
                      mch.y,
                      mch.phi,
                      mch.tgl,
                      mch.signed1Pt,
                      mch.pt,
                      mch.eta,
                      mft.z,
                      mft.x,
                      mft.y,
                      mft.phi,
                      mft.tgl,
                      mft.signed1Pt,
                      mft.pt,
                      mft.eta,
                      globalAtMatchingValues.z,
                      globalAtMatchingValues.x,
                      globalAtMatchingValues.y,
                      globalAtMatchingValues.phi,
                      globalAtMatchingValues.tgl,
                      globalAtMatchingValues.signed1Pt,
                      globalAtMatchingValues.pt,
                      globalAtMatchingValues.eta,
                      mchTruth.pt,
                      mchTruth.eta,
                      mchTruth.phi,
                      mftTruth.pt,
                      mftTruth.eta,
                      mftTruth.phi,
                      globalTruth.pt,
                      globalTruth.eta,
                      globalTruth.phi,
                      deltaX,
                      deltaY,
                      deltaPhi,
                      deltaTgl,
                      deltaSigned1Pt,
                      matchingChi2,
                      mc.isCorrect,
                      mc.category,
                      mc.mchPdgCode,
                      mc.mftPdgCode,
                      mc.mchFirstMotherPdgCode);

        fillHistograms(HistogramGroup::All, mch, mft, *globalKinematics, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        if constexpr (isMC) {
          fillMcCategoryHistograms(mc.category, mch, mft, *globalKinematics, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        }
      }
    }

    if constexpr (isMC) {
      for (const auto& candidateCount : outputCandidatesPerMcMuon) {
        registry.fill(HIST("mc/nOutputCandidatesPerMchMuon"), candidateCount.second);
      }
    }

  }

  void processData(aod::BCsWithTimestamps const& bcs,
                   Collisions const& collisions,
                   FwdTracksData const& fwdTracks,
                   MftTracksData const& mftTracks)
  {
    prepareRun(bcs);
    const std::unordered_map<int64_t, std::unordered_set<int64_t>> emptyMap;
    runStudy<false>(fwdTracks, mftTracks, collisions, emptyMap);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processData, "Process real data", true);

  void processDataWithFTTCA(aod::BCsWithTimestamps const& bcs,
                            Collisions const& collisions,
                            FwdTracksData const& fwdTracks,
                            MftTracksData const& mftTracks,
                            aod::FwdTrackAssoc const& fwdTrackAssociations)
  {
    prepareRun(bcs);
    const std::unordered_map<int64_t, std::unordered_set<int64_t>> emptyMap;
    const auto fwdTrackToCollisions = makeFwdTrackToCollisions(fwdTrackAssociations);
    runStudy<false>(fwdTracks, mftTracks, collisions, emptyMap, &fwdTrackToCollisions);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processDataWithFTTCA, "Process real data with FwdTrackAssoc", false);

  void processMC(aod::BCsWithTimestamps const& bcs,
                 Collisions const& collisions,
                 FwdTracksMC const& fwdTracks,
                 MftTracksMC const& mftTracks,
                 aod::McParticles const& /*mcParticles*/)
  {
    prepareRun(bcs);
    std::unordered_map<int64_t, std::unordered_set<int64_t>> mftParticlesByMcCollision;
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_mcParticle()) {
        continue;
      }
      const auto& particle = mftTrack.mcParticle();
      if (particle.mcCollisionId() >= 0) {
        mftParticlesByMcCollision[particle.mcCollisionId()].insert(particle.globalIndex());
      }
    }
    runStudy<true>(fwdTracks, mftTracks, collisions, mftParticlesByMcCollision);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processMC, "Process MC data", false);

  void processMCWithFTTCA(aod::BCsWithTimestamps const& bcs,
                          Collisions const& collisions,
                          FwdTracksMC const& fwdTracks,
                          MftTracksMC const& mftTracks,
                          aod::McParticles const& /*mcParticles*/,
                          aod::FwdTrackAssoc const& fwdTrackAssociations)
  {
    prepareRun(bcs);
    std::unordered_map<int64_t, std::unordered_set<int64_t>> mftParticlesByMcCollision;
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_mcParticle()) {
        continue;
      }
      const auto& particle = mftTrack.mcParticle();
      if (particle.mcCollisionId() >= 0) {
        mftParticlesByMcCollision[particle.mcCollisionId()].insert(particle.globalIndex());
      }
    }
    const auto fwdTrackToCollisions = makeFwdTrackToCollisions(fwdTrackAssociations);
    runStudy<true>(fwdTracks, mftTracks, collisions, mftParticlesByMcCollision, &fwdTrackToCollisions);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processMCWithFTTCA, "Process MC data with FwdTrackAssoc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalTrackMatchingStudy>(cfgc)};
}
