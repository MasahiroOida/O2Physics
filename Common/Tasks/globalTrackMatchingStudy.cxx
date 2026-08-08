// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file globalTrackMatchingStudy.cxx
/// \brief Matching-plane residual study for reconstructed global muon tracks.

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/Tasks/mftMchMatchingCategories.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
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
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TAxis.h>
#include <TGeoManager.h>
#include <TH1.h>

#include <cmath>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;

namespace matchingcategory = o2::analysis::mftmchmatching;

namespace o2::aod::globaltrackmatchingstudy
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(GlobalTrackId, globalTrackId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
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

  Produces<aod::GlobalTrackMatchingStudies> matchingStudies;

  Configurable<float> matchingPlaneZ{"cfgMatchingPlaneZ", -77.5f, "MFT-MCH matching-plane z (cm)"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of GRPMagField"};
  Configurable<std::string> geometryPath{"geometryPath", "GLO/Config/GeometryAligned", "CCDB path of aligned geometry"};

  ConfigurableAxis axisX{"axisX", {400, -20., 20.}, "#it{x} (cm)"};
  ConfigurableAxis axisY{"axisY", {400, -20., 20.}, "#it{y} (cm)"};
  ConfigurableAxis axisDeltaX{"axisDeltaX", {500, -5., 5.}, "#Delta #it{x} = #it{x}_{MFT} - #it{x}_{MCH} (cm)"};
  ConfigurableAxis axisDeltaY{"axisDeltaY", {500, -5., 5.}, "#Delta #it{y} = #it{y}_{MFT} - #it{y}_{MCH} (cm)"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {500, -0.5, 0.5}, "#Delta #it{#varphi} = #it{#varphi}_{MFT} - #it{#varphi}_{MCH} (rad)"};
  ConfigurableAxis axisDeltaTgl{"axisDeltaTgl", {500, -0.5, 0.5}, "#Delta tan(#it{#lambda})"};
  ConfigurableAxis axisDeltaSigned1Pt{"axisDeltaSigned1Pt", {500, -10., 10.}, "#Delta(#it{q}/#it{p}_{T}) (#it{c}/GeV)"};
  ConfigurableAxis axisMatchingChi2{"axisMatchingChi2", {500, 0., 100.}, "matching #chi^{2}"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::globaltracking::MatchGlobalFwd matching;
  int runNumber = -1;
  float bzAtMftCenter = 0.f;

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
  }

  void init(InitContext const&)
  {
    if (doprocessData == doprocessMC) {
      LOGF(fatal, "Enable exactly one of processData and processMC");
    }
    if (!std::isfinite(matchingPlaneZ.value)) {
      LOGF(fatal, "cfgMatchingPlaneZ must be finite");
    }

    ccdbManager->setURL(ccdbUrl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();

    addStudyHistograms("all/");
    if (doprocessMC) {
      addStudyHistograms("mc/correctPrimaryMuon/");
      addStudyHistograms("mc/correctPiKDecayMuon/");
      addStudyHistograms("mc/correctOtherMuon/");
      addStudyHistograms("mc/falseWithTrueCounterpart/");
      addStudyHistograms("mc/noTrueMftCounterpart/");
      auto categoryHistogram = registry.add<TH1>("mc/category", "Global-track MC matching category;category;tracks", HistType::kTH1F, {{5, -0.5, 4.5}});
      categoryHistogram->GetXaxis()->SetBinLabel(1, "correct primary muon");
      categoryHistogram->GetXaxis()->SetBinLabel(2, "correct #pi/K decay muon");
      categoryHistogram->GetXaxis()->SetBinLabel(3, "correct other muon");
      categoryHistogram->GetXaxis()->SetBinLabel(4, "false with true counterpart");
      categoryHistogram->GetXaxis()->SetBinLabel(5, "no true MFT counterpart");
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
  registry.fill(HIST(prefix "deltaYVsMatchingChi2"), matchingChi2, deltaY)

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
        fillHistograms(HistogramGroup::CorrectPrimaryMuon, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::CorrectPiKDecayMuon:
        fillHistograms(HistogramGroup::CorrectPiKDecayMuon, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::CorrectOtherMuon:
        fillHistograms(HistogramGroup::CorrectOtherMuon, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::FalseWithTrueCounterpart:
        fillHistograms(HistogramGroup::FalseWithTrueCounterpart, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      case Category::NoTrueMftCounterpart:
        fillHistograms(HistogramGroup::NoTrueMftCounterpart, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
        break;
      default:
        break;
    }
  }

  template <bool isMC, typename TFwdTracks, typename TMftTracks>
  void runStudy(const TFwdTracks& fwdTracks,
                const TMftTracks& mftTracks,
                const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftParticlesByMcCollision)
  {
    // matchMCHTrackId stores the global index of the MCH-MID standalone input.
    // Build the lookup explicitly instead of resolving it as an AOD relation: the
    // latter can refer to a non-standalone forward-track row in this input layout.
    std::unordered_map<int64_t, o2::track::TrackParCovFwd> mchTracksByGlobalIndex;
    std::unordered_map<int64_t, McOrigin> mchOriginsByGlobalIndex;
    for (const auto& fwdTrack : fwdTracks) {
      if (fwdTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        continue;
      }
      mchTracksByGlobalIndex.emplace(fwdTrack.globalIndex(), o2::aod::fwdtrackutils::getTrackParCovFwdShift(fwdTrack, 0.f, fwdTrack));
      if constexpr (isMC) {
        mchOriginsByGlobalIndex.emplace(fwdTrack.globalIndex(), mcOrigin(fwdTrack));
      }
    }

    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        continue;
      }

      const int64_t mchTrackId = globalTrack.matchMCHTrackId();
      const int64_t mftTrackId = globalTrack.matchMFTTrackId();
      const auto mchTrack = mchTracksByGlobalIndex.find(mchTrackId);
      if (mchTrack == mchTracksByGlobalIndex.end() || !globalTrack.has_matchMFTTrack()) {
        LOGF(warning, "Skipping global track %lld without a type-3 MCH or valid MFT matching relation: MCH=%lld, MFT=%lld",
             static_cast<long long>(globalTrack.globalIndex()), static_cast<long long>(mchTrackId), static_cast<long long>(mftTrackId));
        continue;
      }

      const auto mftTrack = globalTrack.template matchMFTTrack_as<TMftTracks>();

      // MFT covariance rows are sparse and indexed by matchMFTTrackId. The coordinates
      // and helix propagation used here do not require them, so avoid an invalid row-wise join.
      const auto mftRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(mftTrack, 0.f);
      const auto mchAtMatching = propagateMchToMatchingPlane(mchTrack->second);
      if (!mchAtMatching) {
        LOGF(warning, "MCH propagation failed for global track %lld", static_cast<long long>(globalTrack.globalIndex()));
        continue;
      }
      const auto mftAtMatching = propagateMftToMatchingPlane(mftRaw);
      const auto mch = values(*mchAtMatching);
      const auto mft = values(mftAtMatching);

      const float deltaX = mft.x - mch.x;
      const float deltaY = mft.y - mch.y;
      const float deltaPhi = wrappedDeltaPhi(mft.phi, mch.phi);
      const float deltaTgl = mft.tgl - mch.tgl;
      const float deltaSigned1Pt = mft.signed1Pt - mch.signed1Pt;
      // This is the matching chi2 stored during reconstruction. It is not recalculated
      // when cfgMatchingPlaneZ differs from the plane used by the reconstruction.
      const float matchingChi2 = globalTrack.chi2MatchMCHMFT();

      McClassification mc;
      if constexpr (isMC) {
        mc = classify(mchOriginsByGlobalIndex.at(mchTrackId), mcOrigin(mftTrack), mftParticlesByMcCollision);
      }

      matchingStudies(globalTrack.has_collision() ? globalTrack.collisionId() : -1,
                      globalTrack.globalIndex(),
                      mchTrackId,
                      mftTrackId,
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

      fillHistograms(HistogramGroup::All, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
      if constexpr (isMC) {
        fillMcCategoryHistograms(mc.category, mch, mft, deltaX, deltaY, deltaPhi, deltaTgl, deltaSigned1Pt, matchingChi2);
      }
    }
  }

  void processData(aod::BCsWithTimestamps const& bcs,
                   FwdTracksData const& fwdTracks,
                   MftTracksData const& mftTracks)
  {
    prepareRun(bcs);
    const std::unordered_map<int64_t, std::unordered_set<int64_t>> emptyMap;
    runStudy<false>(fwdTracks, mftTracks, emptyMap);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processData, "Process real data", true);

  void processMC(aod::BCsWithTimestamps const& bcs,
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
    runStudy<true>(fwdTracks, mftTracks, mftParticlesByMcCollision);
  }
  PROCESS_SWITCH(GlobalTrackMatchingStudy, processMC, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalTrackMatchingStudy>(cfgc)};
}
