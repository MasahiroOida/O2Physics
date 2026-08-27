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
//
/// \file mftMchPairProducer.cxx
/// \brief Pair-level MFT-MCH study tables for same-event and mixed-event analyses

#include "PWGDQ/Core/VarManager.h"

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tasks/mftMchMatchingCategories.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TMath.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMFTs = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels, aod::MFTTrkCompColls>;

namespace matchingcategory = o2::analysis::mftmchmatching;

namespace o2::aod
{
namespace mftmchmatchpairs
{
DECLARE_SOA_COLUMN(IsMixedEvent, isMixedEvent, bool);
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, int);
DECLARE_SOA_COLUMN(McMatchCategory, mcMatchCategory, int);
DECLARE_SOA_COLUMN(McMatchStudyCategory, mcMatchStudyCategory, int);
DECLARE_SOA_COLUMN(McPdgCodeMCH, mcPdgCodeMCH, int);
DECLARE_SOA_COLUMN(McPdgCodeMFT, mcPdgCodeMFT, int);
DECLARE_SOA_COLUMN(McFirstMotherPdgCodeMCH, mcFirstMotherPdgCodeMCH, int);
DECLARE_SOA_COLUMN(HasTrueMftForMchMuon, hasTrueMftForMchMuon, bool);
DECLARE_SOA_COLUMN(HasTrueMftForMchMuonInMcCollision, hasTrueMftForMchMuonInMcCollision, bool);
DECLARE_SOA_COLUMN(HasMotherMftForMchMuonInMcCollision, hasMotherMftForMchMuonInMcCollision, bool);
DECLARE_SOA_COLUMN(HasNoTrueOrMotherMftForMchMuonInMcCollision, hasNoTrueOrMotherMftForMchMuonInMcCollision, bool);
DECLARE_SOA_COLUMN(SameMcCollision, sameMcCollision, bool);
DECLARE_SOA_COLUMN(MchMcCollisionId, mchMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
// 1: correct, 0: wrong, -1: truth assignment cannot be evaluated.
DECLARE_SOA_COLUMN(CollisionAssignmentStatus, collisionAssignmentStatus, int);
DECLARE_SOA_COLUMN(NSelectedMchMidMuonsInCollision, nSelectedMchMidMuonsInCollision, int);
DECLARE_SOA_COLUMN(MftMcCollisionId, mftMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MftCollisionId, mftCollisionId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(TimeMCH, timeMCH, double);       //! Absolute MCH time in ns: trackTime + globalBC * bunch spacing
DECLARE_SOA_COLUMN(TimeResMCH, timeResMCH, float); //! MCH track-time resolution in ns
DECLARE_SOA_COLUMN(TimeMFT, timeMFT, double);       //! Absolute MFT time in ns: trackTime + globalBC * bunch spacing
DECLARE_SOA_COLUMN(TimeResMFT, timeResMFT, float); //! MFT track-time resolution in ns
DECLARE_SOA_COLUMN(ZMCH, zMCH, float);
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TglMCH, tglMCH, float);
DECLARE_SOA_COLUMN(InvQPtMCH, invQPtMCH, float);
DECLARE_SOA_COLUMN(PtMCH, ptMCH, float);
DECLARE_SOA_COLUMN(EtaMCH, etaMCH, float);
DECLARE_SOA_COLUMN(PMCH, pMCH, float);
DECLARE_SOA_COLUMN(Chi2MCH, chi2MCH, float);
DECLARE_SOA_COLUMN(InvQPtMCHAtMftLastDisk, invQPtMCHAtMftLastDisk, float);
DECLARE_SOA_COLUMN(PtMCHAtMftLastDisk, ptMCHAtMftLastDisk, float);
DECLARE_SOA_COLUMN(ZMFT, zMFT, float);
DECLARE_SOA_COLUMN(XMFT, xMFT, float);
DECLARE_SOA_COLUMN(YMFT, yMFT, float);
DECLARE_SOA_COLUMN(PhiMFT, phiMFT, float);
DECLARE_SOA_COLUMN(TglMFT, tglMFT, float);
DECLARE_SOA_COLUMN(InvQPtMFT, invQPtMFT, float);
DECLARE_SOA_COLUMN(PtMFT, ptMFT, float);
DECLARE_SOA_COLUMN(EtaMFT, etaMFT, float);
DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);
DECLARE_SOA_COLUMN(NClustersMFT, nClustersMFT, int);
DECLARE_SOA_COLUMN(InvQPtMFTOwn, invQPtMFTOwn, float);
DECLARE_SOA_COLUMN(PtMFTOwn, ptMFTOwn, float);
DECLARE_SOA_COLUMN(EtaMFTOwn, etaMFTOwn, float);
DECLARE_SOA_COLUMN(InvQPtMFTAtMftLastDisk, invQPtMFTAtMftLastDisk, float);
DECLARE_SOA_COLUMN(PtMFTAtMftLastDisk, ptMFTAtMftLastDisk, float);
DECLARE_SOA_COLUMN(PhiMFTPVWithMchMomentum, phiMFTPVWithMchMomentum, float);
DECLARE_SOA_COLUMN(EtaMFTPVWithMchMomentum, etaMFTPVWithMchMomentum, float);
DECLARE_SOA_COLUMN(PMCHAtMftZ, pMCHAtMftZ, float);
DECLARE_SOA_COLUMN(PMCHMinusPPrime, pMCHMinusPPrime, float);
} // namespace mftmchmatchpairs

namespace mchmidvertextracks
{
DECLARE_SOA_COLUMN(ZMCHVertex, zMCHVertex, float);
DECLARE_SOA_COLUMN(XMCHVertex, xMCHVertex, float);
DECLARE_SOA_COLUMN(YMCHVertex, yMCHVertex, float);
DECLARE_SOA_COLUMN(PhiMCHVertex, phiMCHVertex, float);
DECLARE_SOA_COLUMN(TglMCHVertex, tglMCHVertex, float);
DECLARE_SOA_COLUMN(InvQPtMCHVertex, invQPtMCHVertex, float);
DECLARE_SOA_COLUMN(PtMCHVertex, ptMCHVertex, float);
DECLARE_SOA_COLUMN(EtaMCHVertex, etaMCHVertex, float);
DECLARE_SOA_COLUMN(PMCHVertex, pMCHVertex, float);
} // namespace mchmidvertextracks

namespace mftstandalonetracks
{
DECLARE_SOA_COLUMN(TrackTimeMFT, trackTimeMFT, float); //! MFT standalone track time relative to its source collision BC, in ns
DECLARE_SOA_COLUMN(TimeResMFT, timeResMFT, float);     //! MFT standalone track-time resolution, in ns
} // namespace mftstandalonetracks

namespace mchmidstandalonetracks
{
DECLARE_SOA_COLUMN(TrackTimeMCH, trackTimeMCH, float); //! MCH-MID standalone track time relative to its source collision BC, in ns
DECLARE_SOA_COLUMN(TimeResMCH, timeResMCH, float);     //! MCH-MID standalone track-time resolution, in ns
} // namespace mchmidstandalonetracks

DECLARE_SOA_TABLE(MchMidVertexTracks, "AOD", "MCHMIDVTX",
                  o2::soa::Index<>,
                  mftmchmatchpairs::MchCollisionId,
                  mftmchmatchpairs::MchTrackId,
                  mchmidvertextracks::ZMCHVertex,
                  mchmidvertextracks::XMCHVertex,
                  mchmidvertextracks::YMCHVertex,
                  mchmidvertextracks::PhiMCHVertex,
                  mchmidvertextracks::TglMCHVertex,
                  mchmidvertextracks::InvQPtMCHVertex,
                  mchmidvertextracks::PtMCHVertex,
                  mchmidvertextracks::EtaMCHVertex,
                  mchmidvertextracks::PMCHVertex);

DECLARE_SOA_TABLE(MftStandaloneTracks, "AOD", "MFTSTAND",
                  o2::soa::Index<>,
                  mftmchmatchpairs::MftMcCollisionId,
                  mftmchmatchpairs::MftCollisionId,
                  mftmchmatchpairs::CollisionAssignmentStatus,
                  mftmchmatchpairs::MftTrackId,
                  mftstandalonetracks::TrackTimeMFT,
                  mftstandalonetracks::TimeResMFT);

DECLARE_SOA_TABLE(MchMidStandaloneTracks, "AOD", "MCHMIDST",
                  o2::soa::Index<>,
                  mftmchmatchpairs::MchMcCollisionId,
                  mftmchmatchpairs::MchCollisionId,
                  mftmchmatchpairs::CollisionAssignmentStatus,
                  mftmchmatchpairs::MchTrackId,
                  mchmidstandalonetracks::TrackTimeMCH,
                  mchmidstandalonetracks::TimeResMCH);

namespace mftmchmatchpairs
{
DECLARE_SOA_INDEX_COLUMN(MchMidVertexTrack, mchMidVertexTrack); //! Associated MCHMIDVTX row
DECLARE_SOA_INDEX_COLUMN(MftStandaloneTrack, mftStandaloneTrack); //! Associated MFTSTAND row
DECLARE_SOA_INDEX_COLUMN(MchMidStandaloneTrack, mchMidStandaloneTrack); //! Associated MCHMIDST row
} // namespace mftmchmatchpairs

DECLARE_SOA_TABLE(MftMchMatchPairs, "AOD", "MFTMCHPAIR",
                  mftmchmatchpairs::IsMixedEvent,
                  mftmchmatchpairs::IsCorrectMatch,
                  mftmchmatchpairs::McMatchCategory,
                  mftmchmatchpairs::McMatchStudyCategory,
                  mftmchmatchpairs::McPdgCodeMCH,
                  mftmchmatchpairs::McPdgCodeMFT,
                  mftmchmatchpairs::McFirstMotherPdgCodeMCH,
                  mftmchmatchpairs::HasTrueMftForMchMuon,
                  mftmchmatchpairs::HasTrueMftForMchMuonInMcCollision,
                  mftmchmatchpairs::HasMotherMftForMchMuonInMcCollision,
                  mftmchmatchpairs::HasNoTrueOrMotherMftForMchMuonInMcCollision,
                  mftmchmatchpairs::SameMcCollision,
                  mftmchmatchpairs::MchMcCollisionId,
                  mftmchmatchpairs::MchCollisionId,
                  mftmchmatchpairs::NSelectedMchMidMuonsInCollision,
                  mftmchmatchpairs::MftMcCollisionId,
                  mftmchmatchpairs::MftCollisionId,
                  mftmchmatchpairs::MchTrackId,
                  mftmchmatchpairs::MftTrackId,
                  mftmchmatchpairs::MchMidVertexTrackId,
                  mftmchmatchpairs::MftStandaloneTrackId,
                  mftmchmatchpairs::MchMidStandaloneTrackId,
                  mftmchmatchpairs::TimeMCH,
                  mftmchmatchpairs::TimeResMCH,
                  mftmchmatchpairs::TimeMFT,
                  mftmchmatchpairs::TimeResMFT,
                  mftmchmatchpairs::ZMCH,
                  mftmchmatchpairs::XMCH,
                  mftmchmatchpairs::YMCH,
                  mftmchmatchpairs::PhiMCH,
                  mftmchmatchpairs::TglMCH,
                  mftmchmatchpairs::InvQPtMCH,
                  mftmchmatchpairs::PtMCH,
                  mftmchmatchpairs::EtaMCH,
                  mftmchmatchpairs::PMCH,
                  mftmchmatchpairs::Chi2MCH,
                  mftmchmatchpairs::InvQPtMCHAtMftLastDisk,
                  mftmchmatchpairs::PtMCHAtMftLastDisk,
                  mftmchmatchpairs::ZMFT,
                  mftmchmatchpairs::XMFT,
                  mftmchmatchpairs::YMFT,
                  mftmchmatchpairs::PhiMFT,
                  mftmchmatchpairs::TglMFT,
                  mftmchmatchpairs::InvQPtMFT,
                  mftmchmatchpairs::PtMFT,
                  mftmchmatchpairs::EtaMFT,
                  mftmchmatchpairs::Chi2MFT,
                  mftmchmatchpairs::NClustersMFT,
                  mftmchmatchpairs::InvQPtMFTOwn,
                  mftmchmatchpairs::PtMFTOwn,
                  mftmchmatchpairs::EtaMFTOwn,
                  mftmchmatchpairs::InvQPtMFTAtMftLastDisk,
                  mftmchmatchpairs::PtMFTAtMftLastDisk,
                  mftmchmatchpairs::PhiMFTPVWithMchMomentum,
                  mftmchmatchpairs::EtaMFTPVWithMchMomentum,
                  mftmchmatchpairs::PMCHAtMftZ,
                  mftmchmatchpairs::PMCHMinusPPrime);
} // namespace o2::aod

struct mftMchPairProducer {
  Produces<o2::aod::MftMchMatchPairs> mftMchMatchPairs;
  Produces<o2::aod::MchMidVertexTracks> mchMidVertexTracks;
  Produces<o2::aod::MftStandaloneTracks> mftStandaloneTracks;
  Produces<o2::aod::MchMidStandaloneTracks> mchMidStandaloneTracks;
  ////   Variables for selecting muon tracks
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> fPhiMchLow{"cfgPhiMchLow", -TMath::Pi(), ""};
  Configurable<float> fPhiMchUp{"cfgPhiMchUp", TMath::Pi(), ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMFTLow{"cfgEtaMFTlow", -3.6f, ""};
  Configurable<float> fEtaMFTUp{"cfgEtaMFTup", -2.5f, ""};
  Configurable<float> fPhiMFTLow{"cfgPhiMFTLow", -TMath::Pi(), ""};
  Configurable<float> fPhiMFTUp{"cfgPhiMFTUp", TMath::Pi(), ""};
  Configurable<float> fPtMFTLow{"cfgPtMFTLow", 0.1f, ""};

  ////   Variables for selecting collisions
  Configurable<bool> fApplySel8{"cfgApplySel8", false, "Apply sel8 event selection to same-event and mixed-event MCH-MID to MFT histograms"};
  Configurable<bool> fApplyVtxZCut{"cfgApplyVtxZCut", false, "Apply vertex-z event selection to same-event and mixed-event MCH-MID to MFT histograms"};
  Configurable<float> fVtxZMax{"cfgVtxZMax", 10.f, "Maximum absolute collision vertex z"};

  ////   Variables for matching configuration
  Configurable<float> fzMatching{"cfgzMatching", -77.5f, "Plane for MFT-MCH matching"};
  Configurable<bool> fFillAllMchMidMftHelixDeltas{"cfgFillAllMchMidMftHelixDeltas", true, "Process all same-collision MCH-MID/MFT pairs at the matching plane for MFTMCHPAIR output"};
  Configurable<bool> fFillMixedMchMidMftHelixDeltas{"cfgFillMixedMchMidMftHelixDeltas", true, "Process mixed-event MCH-MID/MFT pairs at the matching plane for MFTMCHPAIR output"};
  Configurable<bool> fFillQaHistograms{"cfgFillQaHistograms", false, "Fill QA/debug histograms not needed by the pair-table based downstream analysis"};
  Configurable<bool> fUseBestTimeCollision{"cfgUseBestTimeCollision", false, "Use the time-compatible collision with the smallest normalized time residual for MCH-MID tracks"};
  Configurable<bool> fUseMcTruthCollisionAssociation{"cfgUseMcTruthCollisionAssociation", true,
                                                     "For MC only, require a unique selected reconstructed collision with the same MC collision label as the track"};
  Configurable<bool> fRequireSingleMchCompatibleCollision{"cfgRequireSingleMchCompatibleCollision", false, "Analyze only MCH-MID tracks with exactly one time-compatible collision"};
  Configurable<int> fEventMixingDepth{"cfgEventMixingDepth", 5, "Number of collisions per vertex-z bin used for event mixing"};
  Configurable<int64_t> fMinMixingBcGap{"cfgMinMixingBcGap", 0, "Minimum absolute global-BC gap between collisions used for MCH-MID/MFT event mixing"};
  ConfigurableAxis fEventMixingVtxBins{"cfgEventMixingVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Vertex-z bins for event mixing"};

  Configurable<float> fSamplingFraction{"cfgSamplingFraction", 1.f, "Fraction of randomly selected events to be processed"};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  double mBzAtMftCenter{0};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field
  Service<o2::ccdb::BasicCCDBManager> ccdbManager{};

  o2::parameters::GRPMagField* fGrpMag = nullptr;

  std::uniform_real_distribution<double> mDistribution{0.0, 1.0};
  std::mt19937 mGenerator;

  o2::globaltracking::MatchGlobalFwd mMatching;

  std::unordered_map<int64_t, int32_t> mftCovIndexes;

  HistogramRegistry registry{"registry", {}};
  
  struct PropagatedTrackAtMatching {
    float z = 0.f;
    float x = 0.f;
    float y = 0.f;
    float phi = 0.f;
    float tgl = 0.f;
    float invQPt = 0.f;
    float pt = 0.f;
    float eta = 0.f;
  };

  struct PropagatedMchMidTrackAtMatching {
    o2::track::TrackParCovFwd rawTrack;
    PropagatedTrackAtMatching values;
    int64_t globalIndex = -1;
    int64_t mchMidVertexTrackId = -1;
    int64_t mchMidStandaloneTrackId = -1;
    double time = 0.;
    float timeRes = 0.f;
    float chi2 = 0.f;
  };

  struct McTrackOriginInfo {
    bool hasMcParticle = false;
    int64_t mcParticleId = -1;
    int64_t mcCollisionId = -1;
    int pdgCode = 0;
    bool isMuon = false;
    bool isPhysicalPrimary = false;
    bool hasPiKMother = false;
    std::unordered_set<int64_t> piKMotherMcParticleIds;
    bool hasMotherMcParticle = false;
    int64_t motherMcParticleId = -1;
    int firstMotherPdgCode = 0;
  };

  struct MftTrackForPairing {
    o2::track::TrackParCovFwd track;
    // The own-momentum parameters at the MFT last disk are independent of the MCH track, so they are cached once per MFT track.
    PropagatedTrackAtMatching valuesAtMftLastDisk;
    int64_t globalIndex = -1;
    int64_t mftStandaloneTrackId = -1;
    double time = 0.;
    float timeRes = 0.f;
    float chi2 = 0.f;
    int nClusters = 0;
    bool hasMcParticle = false;
    int64_t mcParticleId = -1;
    McTrackOriginInfo mcOrigin;
  };

  struct PairPropagationResult {
    PropagatedTrackAtMatching mchValuesAtMftLastDisk;
    PropagatedTrackAtMatching mftValuesAtMatching;
    float mchMomentumAtMftZ = 0.f;
    float mchMomentumLoss = 0.f;
  };

  struct McPairInfo {
    int isCorrectMatch = -1;
    int matchCategory = matchingcategory::serialized(matchingcategory::McMatchCategory::Unknown);
    int matchStudyCategory = matchingcategory::serialized(matchingcategory::McMatchStudyCategory::Unknown);
    int mchMcPdgCode = 0;
    int mftMcPdgCode = 0;
    int mchFirstMotherPdgCode = 0;
    bool hasTrueMftForMchMuon = false;
    bool hasTrueMftForMchMuonInMcCollision = false;
    bool hasMotherMftForMchMuonInMcCollision = false;
    bool hasNoTrueOrMotherMftForMchMuonInMcCollision = false;
    bool sameMcCollision = false;
    int64_t mchMcCollisionId = -1;
    int64_t mftMcCollisionId = -1;
  };

  struct BestTimeCollisionChoice {
    int64_t selectedCollisionId = -1;
    int64_t originalCollisionId = -1;
    int nCompatibleCollisions = 0;
    float bestTimeResidual = std::numeric_limits<float>::quiet_NaN();
    float secondBestTimeResidual = std::numeric_limits<float>::quiet_NaN();
  };

  static float etaFromTgl(float tgl)
  {
    return std::asinh(tgl);
  }

  static float ptFromInvQPt(float invQPt)
  {
    constexpr float minAbsInvQPt = 1.e-12f;
    return std::abs(invQPt) > minAbsInvQPt ? std::abs(1.f / invQPt) : 0.f;
  }

  static float pFromPtTgl(float pt, float tgl)
  {
    return pt * std::sqrt(1.f + tgl * tgl);
  }

  template <class TBCS, class TTRACK, class TCOLLS>
  static double absoluteTrackTime(const TTRACK& track, const TCOLLS& collisions)
  {
    if (!track.has_collision()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const auto& collision = collisions.rawIteratorAt(track.collisionId());
    const auto& bc = collision.template bc_as<TBCS>();
    return static_cast<double>(track.trackTime()) +
           static_cast<double>(bc.globalBC()) * o2::constants::lhc::LHCBunchSpacingNS;
  }

  static PropagatedTrackAtMatching makePropagatedTrackValues(const o2::track::TrackParCovFwd& track)
  {
    PropagatedTrackAtMatching values;
    values.z = track.getZ();
    values.x = track.getX();
    values.y = track.getY();
    values.phi = track.getPhi();
    values.tgl = track.getTgl();
    values.invQPt = track.getInvQPt();
    values.pt = ptFromInvQPt(values.invQPt);
    values.eta = etaFromTgl(values.tgl);
    return values;
  }

  std::optional<o2::track::TrackParCovFwd> propagateMchToZWithEnergyLoss(const o2::track::TrackParCovFwd& fwdTrack, float z)
  {
    constexpr float absorberBeginZ = -90.f;
    o2::dataformats::GlobalFwdTrack globalTrack;
    globalTrack.setParameters(fwdTrack.getParameters());
    globalTrack.setZ(fwdTrack.getZ());
    globalTrack.setCovariances(fwdTrack.getCovariances());

    auto mchTrack = mMatching.FwdtoMCH(globalTrack);
    // For planes downstream of the absorber entrance, propagate through the absorber with energy-loss handling
    // but without the Branson correction. Upstream MCH-side planes use the standard covariance extrapolation.
    const bool propagationOk = z >= absorberBeginZ ? o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z)
                                                   : o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    if (!propagationOk) {
      return std::nullopt;
    }
    o2::track::TrackParCovFwd propagatedTrack = mMatching.MCHtoFwd(mchTrack);
    return propagatedTrack;
  }

  template <class C>
  std::optional<o2::track::TrackParCovFwd> propagateMchToCollisionVertex(const o2::track::TrackParCovFwd& fwdTrack, const C& collision)
  {
    o2::dataformats::GlobalFwdTrack globalTrack;
    globalTrack.setParameters(fwdTrack.getParameters());
    globalTrack.setZ(fwdTrack.getZ());
    globalTrack.setCovariances(fwdTrack.getCovariances());

    auto mchTrack = mMatching.FwdtoMCH(globalTrack);
    const bool propagationOk = o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    if (!propagationOk) {
      return std::nullopt;
    }
    o2::track::TrackParCovFwd propagatedTrack = mMatching.MCHtoFwd(mchTrack);
    return propagatedTrack;
  }

  static bool setTrackTotalMomentum(o2::track::TrackParCovFwd& track, float totalMomentum)
  {
    if (!std::isfinite(totalMomentum) || totalMomentum <= 0.f) {
      return false;
    }
    const auto tgl = static_cast<float>(track.getTgl());
    const auto pt = totalMomentum / std::sqrt(1.f + tgl * tgl);
    if (!std::isfinite(pt) || pt <= 0.f) {
      return false;
    }
    const auto chargeSign = std::signbit(track.getInvQPt()) ? -1.f : 1.f;
    track.setInvQPt(chargeSign / pt);
    return true;
  }

  std::optional<PropagatedTrackAtMatching> propagateMftToMatchingWithMchMomentum(const o2::track::TrackParCovFwd& mftTrack,
                                                                                 float mchMomentumAtMftZ) const
  {
    // Pair-dependent MFT propagation: replace only the momentum magnitude with the MCH momentum propagated to the MFT z.
    auto mftWithMchMomentum = mftTrack;
    if (!setTrackTotalMomentum(mftWithMchMomentum, mchMomentumAtMftZ)) {
      return std::nullopt;
    }
    mftWithMchMomentum.propagateToZhelix(fzMatching.value, mBzAtMftCenter);
    return makePropagatedTrackValues(mftWithMchMomentum);
  }

  std::optional<PropagatedTrackAtMatching> propagateMftToZWithMchMomentum(const o2::track::TrackParCovFwd& mftTrack,
                                                                          float mchMomentumAtMftZ,
                                                                          float z) const
  {
    auto mftWithMchMomentum = mftTrack;
    if (!setTrackTotalMomentum(mftWithMchMomentum, mchMomentumAtMftZ)) {
      return std::nullopt;
    }
    mftWithMchMomentum.propagateToZhelix(z, mBzAtMftCenter);
    return makePropagatedTrackValues(mftWithMchMomentum);
  }

  static bool inRange(float value, float low, float up)
  {
    return value >= low && value <= up;
  }

  static bool inPhiRange(float phi, float low, float up)
  {
    if (low <= up) {
      return phi >= low && phi <= up;
    }
    return phi >= low || phi <= up;
  }

  template <class TCOLL>
  bool isSelectedCollision(const TCOLL& collision) const
  {
    if (fApplySel8.value && !collision.sel8()) {
      return false;
    }
    if (fApplyVtxZCut.value && std::abs(collision.posZ()) > fVtxZMax.value) {
      return false;
    }
    return true;
  }

  template <class TTRACK>
  bool hasRequiredMchCompatibleCollisionMultiplicity(const TTRACK& track) const
  {
    return !fRequireSingleMchCompatibleCollision.value || track.compatibleCollIds().size() == 1;
  }

  template <class TBCS, class TTRACK, class TCOLLS>
  BestTimeCollisionChoice chooseBestTimeCollision(const TTRACK& track,
                                                  const TCOLLS& collisions,
                                                  const std::unordered_set<int64_t>& selectedCollisionIds) const
  {
    BestTimeCollisionChoice choice;
    if (!track.has_collision()) {
      return choice;
    }

    choice.originalCollisionId = static_cast<int64_t>(track.collisionId());
    choice.selectedCollisionId = choice.originalCollisionId;

    if (!fUseBestTimeCollision.value) {
      choice.nCompatibleCollisions = 1;
      return choice;
    }

    const auto& compatibleCollisionIds = track.compatibleCollIds();
    choice.nCompatibleCollisions = static_cast<int>(compatibleCollisionIds.size());
    if (compatibleCollisionIds.empty()) {
      return choice;
    }

    const auto& originalCollision = collisions.rawIteratorAt(choice.originalCollisionId);
    const auto& originalBC = originalCollision.template bc_as<TBCS>();
    const int64_t trackBC = static_cast<int64_t>(originalBC.globalBC());
    const float trackTime = track.trackTime();
    const float trackTimeRes = track.trackTimeRes();
    float bestAbsResidual = std::numeric_limits<float>::max();
    float secondBestAbsResidual = std::numeric_limits<float>::max();

    for (const auto compatibleCollisionId : compatibleCollisionIds) {
      const int64_t collisionId = static_cast<int64_t>(compatibleCollisionId);
      if (selectedCollisionIds.find(collisionId) == selectedCollisionIds.end()) {
        continue;
      }

      const auto& candidateCollision = collisions.rawIteratorAt(collisionId);
      const auto& candidateBC = candidateCollision.template bc_as<TBCS>();
      const int64_t collisionBC = static_cast<int64_t>(candidateBC.globalBC());
      const float deltaTime = trackTime - candidateCollision.collisionTime() + (trackBC - collisionBC) * o2::constants::lhc::LHCBunchSpacingNS;
      const float sigmaTime = std::sqrt(trackTimeRes * trackTimeRes + candidateCollision.collisionTimeRes() * candidateCollision.collisionTimeRes());
      const float residual = sigmaTime > 0.f ? deltaTime / sigmaTime : deltaTime;
      const float absResidual = std::abs(residual);

      if (absResidual < bestAbsResidual) {
        secondBestAbsResidual = bestAbsResidual;
        choice.secondBestTimeResidual = choice.bestTimeResidual;
        bestAbsResidual = absResidual;
        choice.bestTimeResidual = residual;
        choice.selectedCollisionId = collisionId;
      } else if (absResidual < secondBestAbsResidual) {
        secondBestAbsResidual = absResidual;
        choice.secondBestTimeResidual = residual;
      }
    }

    return choice;
  }

  template <class TTRACK, class TCOLLS>
  std::optional<int64_t> chooseMcTruthCollision(const TTRACK& track,
                                                 const TCOLLS& collisions,
                                                 const std::unordered_set<int64_t>& selectedCollisionIds) const
  {
    // Follow the PWGDQ MC-association QA convention: a fake/noise label is not
    // usable as a truth collision, even when a MC-particle index is present.
    constexpr uint8_t noiseBit = 0x1 << 6;
    constexpr uint8_t fakeBit = 0x1 << 7;
    if (!track.has_mcParticle() || (track.mcMask() & (noiseBit | fakeBit))) {
      return std::nullopt;
    }

    const int64_t truthMcCollisionId = track.mcParticle().mcCollisionId();
    if (truthMcCollisionId < 0) {
      return std::nullopt;
    }

    int64_t matchingCollisionId = -1;
    int nMatchingCollisions = 0;
    for (const auto& collision : collisions) {
      const int64_t collisionId = static_cast<int64_t>(collision.globalIndex());
      if (selectedCollisionIds.find(collisionId) == selectedCollisionIds.end() || !collision.has_mcCollision() ||
          collision.mcCollisionId() != truthMcCollisionId) {
        continue;
      }
      matchingCollisionId = collisionId;
      ++nMatchingCollisions;
    }

    // A split MC collision has several reconstructed vertices. Do not choose one
    // arbitrarily: this is a truth-association closure selection, not a new
    // reconstruction-level association.
    return nMatchingCollisions == 1 ? std::optional<int64_t>{matchingCollisionId} : std::nullopt;
  }

  template <bool isMc, class TBCS, class TTRACK, class TCOLLS>
  std::optional<BestTimeCollisionChoice> choosePairCollision(const TTRACK& track,
                                                              const TCOLLS& collisions,
                                                              const std::unordered_set<int64_t>& selectedCollisionIds) const
  {
    if constexpr (isMc) {
      if (fUseMcTruthCollisionAssociation.value) {
        const auto truthCollisionId = chooseMcTruthCollision(track, collisions, selectedCollisionIds);
        if (!truthCollisionId) {
          return std::nullopt;
        }
        BestTimeCollisionChoice choice;
        choice.selectedCollisionId = *truthCollisionId;
        choice.originalCollisionId = track.has_collision() ? static_cast<int64_t>(track.collisionId()) : -1;
        choice.nCompatibleCollisions = 1;
        return choice;
      }
    }

    if (!track.has_collision()) {
      return std::nullopt;
    }
    return chooseBestTimeCollision<TBCS>(track, collisions, selectedCollisionIds);
  }

  bool isSelectedMchMidAtMatching(const PropagatedTrackAtMatching& track, float rAbs) const
  {
    if (!inRange(track.eta, fEtaMchLow.value, fEtaMchUp.value)) {
      return false;
    }
    if (!inPhiRange(track.phi, fPhiMchLow.value, fPhiMchUp.value)) {
      return false;
    }
    if (track.pt < fPtMchLow.value) {
      return false;
    }
    if (!inRange(rAbs, fRabsLow.value, fRabsUp.value)) {
      return false;
    }
    return true;
  }

  bool isSelectedMFTAtMatching(const PropagatedTrackAtMatching& track) const
  {
    if (!inRange(track.eta, fEtaMFTLow.value, fEtaMFTUp.value)) {
      return false;
    }
    if (!inPhiRange(track.phi, fPhiMFTLow.value, fPhiMFTUp.value)) {
      return false;
    }
    if (track.pt < fPtMFTLow.value) {
      return false;
    }
    return true;
  }

  std::optional<PairPropagationResult> propagatePair(const o2::track::TrackParCovFwd& mchTrack,
                                                      const MftTrackForPairing& mftTrack,
                                                      bool isMixedEvent)
  {
    auto mchAtMftZ = propagateMchToZWithEnergyLoss(mchTrack, mftTrack.track.getZ());
    if (!mchAtMftZ) {
      return std::nullopt;
    }

    const float mchMomentumAtMftZ = static_cast<float>(mchAtMftZ->getP());
    const float mchMomentumLoss = static_cast<float>(mchTrack.getP() - mchAtMftZ->getP());
    if (fFillQaHistograms) {
      if (isMixedEvent) {
        registry.fill(HIST("mixedMchMidMftHelix/debugPMinusPPrime"), mchMomentumLoss);
      } else {
        registry.fill(HIST("allMchMidMftHelix/debugPMinusPPrime"), mchMomentumLoss);
      }
    }

    auto mftAtMatching = propagateMftToMatchingWithMchMomentum(mftTrack.track, mchMomentumAtMftZ);
    if (!mftAtMatching) {
      return std::nullopt;
    }
    // Keep the existing QA scope: this histogram is filled only while processing mixed-event pairs.
    if (isMixedEvent && fFillQaHistograms) {
      registry.fill(HIST("allMchMidMftHelix/zMftAtMatching"), mftAtMatching->z);
    }
    if (!isSelectedMFTAtMatching(*mftAtMatching)) {
      return std::nullopt;
    }

    return PairPropagationResult{makePropagatedTrackValues(*mchAtMftZ),
                                 *mftAtMatching,
                                 mchMomentumAtMftZ,
                                 mchMomentumLoss};
  }

  template <class TTRACK>
  McTrackOriginInfo getMcTrackOriginInfo(const TTRACK& track) const
  {
    McTrackOriginInfo info;
    if (!track.has_mcParticle()) {
      return info;
    }

    const auto& mcParticle = track.mcParticle();
    info.hasMcParticle = true;
    info.mcParticleId = mcParticle.globalIndex();
    info.mcCollisionId = mcParticle.mcCollisionId();
    info.pdgCode = mcParticle.pdgCode();
    info.isMuon = std::abs(mcParticle.pdgCode()) == 13;
    info.isPhysicalPrimary = mcParticle.isPhysicalPrimary();

    const auto& motherParticles = mcParticle.template mothers_as<aod::McParticles>();
    if (!motherParticles.empty()) {
      const auto& motherParticle = motherParticles[0];
      info.hasMotherMcParticle = true;
      info.motherMcParticleId = motherParticle.globalIndex();
      info.firstMotherPdgCode = motherParticle.pdgCode();
    }
    for (const auto& motherParticle : motherParticles) {
      const int absMotherPdg = std::abs(motherParticle.pdgCode());
      if (absMotherPdg == 211 || absMotherPdg == 321) {
        info.hasPiKMother = true;
        info.piKMotherMcParticleIds.insert(motherParticle.globalIndex());
      }
    }
    return info;
  }

  bool hasMftMcParticleInMcCollision(const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftMcParticleIdsPerMcCollision,
                                     int64_t mcCollisionId,
                                     int64_t mcParticleId) const
  {
    if (mcCollisionId < 0 || mcParticleId < 0) {
      return false;
    }
    const auto mftMcIdsInMcCollisionIt = mftMcParticleIdsPerMcCollision.find(mcCollisionId);
    return mftMcIdsInMcCollisionIt != mftMcParticleIdsPerMcCollision.end() &&
           mftMcIdsInMcCollisionIt->second.find(mcParticleId) != mftMcIdsInMcCollisionIt->second.end();
  }

  bool hasAnyMftMcParticleInMcCollision(const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftMcParticleIdsPerMcCollision,
                                        int64_t mcCollisionId,
                                        const std::unordered_set<int64_t>& mcParticleIds) const
  {
    if (mcCollisionId < 0 || mcParticleIds.empty()) {
      return false;
    }
    const auto mftMcIdsInMcCollisionIt = mftMcParticleIdsPerMcCollision.find(mcCollisionId);
    if (mftMcIdsInMcCollisionIt == mftMcParticleIdsPerMcCollision.end()) {
      return false;
    }
    for (const auto mcParticleId : mcParticleIds) {
      if (mftMcIdsInMcCollisionIt->second.find(mcParticleId) != mftMcIdsInMcCollisionIt->second.end()) {
        return true;
      }
    }
    return false;
  }

  McPairInfo makeMcPairInfo(const McTrackOriginInfo& mchOrigin,
                            const MftTrackForPairing& mftTrack,
                            int64_t collisionId,
                            const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftMcParticleIdsPerCollision,
                            const std::unordered_map<int64_t, std::unordered_set<int64_t>>& mftMcParticleIdsPerMcCollision) const
  {
    McPairInfo info;
    info.mchMcPdgCode = mchOrigin.pdgCode;
    info.mftMcPdgCode = mftTrack.mcOrigin.pdgCode;
    info.mchFirstMotherPdgCode = mchOrigin.firstMotherPdgCode;
    info.mchMcCollisionId = mchOrigin.hasMcParticle ? mchOrigin.mcCollisionId : -1;
    info.mftMcCollisionId = mftTrack.mcOrigin.hasMcParticle ? mftTrack.mcOrigin.mcCollisionId : -1;
    const bool isSameMcCollision = mchOrigin.hasMcParticle && mftTrack.hasMcParticle && mchOrigin.mcCollisionId >= 0 &&
                                   mftTrack.mcOrigin.mcCollisionId >= 0 && mchOrigin.mcCollisionId == mftTrack.mcOrigin.mcCollisionId;
    const bool isSameMcParticle = mchOrigin.hasMcParticle && mftTrack.hasMcParticle && mchOrigin.mcParticleId == mftTrack.mcParticleId;
    const bool isCandidatePiKMotherMatch = mchOrigin.isMuon && mchOrigin.hasPiKMother && isSameMcCollision &&
                                           mchOrigin.piKMotherMcParticleIds.find(mftTrack.mcParticleId) != mchOrigin.piKMotherMcParticleIds.end();
    const bool isCorrectMatch = matchingcategory::isCorrectPair(isSameMcParticle, isCandidatePiKMotherMatch);

    info.isCorrectMatch = isCorrectMatch ? 1 : 0;
    const auto mftMcIdsIt = mftMcParticleIdsPerCollision.find(collisionId);
    if (mftMcIdsIt != mftMcParticleIdsPerCollision.end() && mchOrigin.hasMcParticle) {
      info.hasTrueMftForMchMuon = mftMcIdsIt->second.find(mchOrigin.mcParticleId) != mftMcIdsIt->second.end();
    }
    if (mchOrigin.hasMcParticle) {
      info.hasTrueMftForMchMuonInMcCollision = hasMftMcParticleInMcCollision(mftMcParticleIdsPerMcCollision,
                                                                            mchOrigin.mcCollisionId,
                                                                            mchOrigin.mcParticleId);
    }
    if (mchOrigin.isMuon && mchOrigin.hasPiKMother) {
      info.hasMotherMftForMchMuonInMcCollision = hasAnyMftMcParticleInMcCollision(mftMcParticleIdsPerMcCollision,
                                                                                   mchOrigin.mcCollisionId,
                                                                                   mchOrigin.piKMotherMcParticleIds);
    }
    info.hasNoTrueOrMotherMftForMchMuonInMcCollision = !info.hasTrueMftForMchMuonInMcCollision && !info.hasMotherMftForMchMuonInMcCollision;
    info.sameMcCollision = isSameMcCollision;
    info.matchCategory = matchingcategory::serialized(
      matchingcategory::classifyPair(isCorrectMatch,
                                     mchOrigin.hasMcParticle,
                                     mftTrack.mcOrigin.hasMcParticle,
                                     mchOrigin.mcCollisionId,
                                     mftTrack.mcOrigin.mcCollisionId));
    info.matchStudyCategory = matchingcategory::serialized(
      matchingcategory::classifyStudyPair(isCorrectMatch,
                                          mchOrigin.isMuon,
                                          mchOrigin.hasPiKMother,
                                          mchOrigin.isPhysicalPrimary,
                                          info.hasTrueMftForMchMuonInMcCollision || info.hasMotherMftForMchMuonInMcCollision));
    return info;
  }

  void fillMftMchPairTable(bool isMixedEvent,
                           const McPairInfo& mcPairInfo,
                           int64_t mchCollisionId,
                           int nSelectedMchMidMuonsInCollision,
                           int64_t mftCollisionId,
                           int64_t mchTrackId,
                           int64_t mftTrackId,
                           int64_t mchMidVertexTrackId,
                           int64_t mftStandaloneTrackId,
                           int64_t mchMidStandaloneTrackId,
                           double timeMCH,
                           float timeResMCH,
                           double timeMFT,
                           float timeResMFT,
                           const PropagatedTrackAtMatching& muonValues,
                           const PropagatedTrackAtMatching& mchValuesAtMftLastDisk,
                           float chi2MCH,
                           const PropagatedTrackAtMatching& mftValues,
                           const PropagatedTrackAtMatching& mftOwnMomentumValues,
                           float chi2MFT,
                           int nClustersMFT,
                           float phiMftPVWithMchMomentum,
                           float etaMftPVWithMchMomentum,
                           float pMchAtMftZ,
                           float pMchMinusPPrime)
  {
    mftMchMatchPairs(
      isMixedEvent,
      mcPairInfo.isCorrectMatch,
      mcPairInfo.matchCategory,
      mcPairInfo.matchStudyCategory,
      mcPairInfo.mchMcPdgCode,
      mcPairInfo.mftMcPdgCode,
      mcPairInfo.mchFirstMotherPdgCode,
      mcPairInfo.hasTrueMftForMchMuon,
      mcPairInfo.hasTrueMftForMchMuonInMcCollision,
      mcPairInfo.hasMotherMftForMchMuonInMcCollision,
      mcPairInfo.hasNoTrueOrMotherMftForMchMuonInMcCollision,
      mcPairInfo.sameMcCollision,
      mcPairInfo.mchMcCollisionId,
      mchCollisionId,
      nSelectedMchMidMuonsInCollision,
      mcPairInfo.mftMcCollisionId,
      mftCollisionId,
      mchTrackId,
      mftTrackId,
      mchMidVertexTrackId,
      mftStandaloneTrackId,
      mchMidStandaloneTrackId,
      timeMCH,
      timeResMCH,
      timeMFT,
      timeResMFT,
      muonValues.z,
      muonValues.x,
      muonValues.y,
      muonValues.phi,
      muonValues.tgl,
      muonValues.invQPt,
      muonValues.pt,
      muonValues.eta,
      pFromPtTgl(muonValues.pt, muonValues.tgl),
      chi2MCH,
      mchValuesAtMftLastDisk.invQPt,
      mchValuesAtMftLastDisk.pt,
      mftValues.z,
      mftValues.x,
      mftValues.y,
      mftValues.phi,
      mftValues.tgl,
      mftValues.invQPt,
      mftValues.pt,
      mftValues.eta,
      chi2MFT,
      nClustersMFT,
      mftOwnMomentumValues.invQPt,
      mftOwnMomentumValues.pt,
      mftOwnMomentumValues.eta,
      mftOwnMomentumValues.invQPt,
      mftOwnMomentumValues.pt,
      phiMftPVWithMchMomentum,
      etaMftPVWithMchMomentum,
      pMchAtMftZ,
      pMchMinusPPrime);
  }

  int64_t fillMchMidVertexTable(int64_t mchCollisionId,
                                int64_t mchTrackId,
                                const PropagatedTrackAtMatching& vertexValues)
  {
    mchMidVertexTracks(
      mchCollisionId,
      mchTrackId,
      vertexValues.z,
      vertexValues.x,
      vertexValues.y,
      vertexValues.phi,
      vertexValues.tgl,
      vertexValues.invQPt,
      vertexValues.pt,
      vertexValues.eta,
      pFromPtTgl(vertexValues.pt, vertexValues.tgl));
    return static_cast<int64_t>(mchMidVertexTracks.lastIndex());
  }

  int64_t fillMftStandaloneTrackTable(int64_t mftCollisionId,
                                      int64_t mftMcCollisionId,
                                      int collisionAssignmentStatus,
                                      int64_t mftTrackId,
                                      float trackTime,
                                      float timeRes)
  {
    mftStandaloneTracks(mftMcCollisionId, mftCollisionId, collisionAssignmentStatus, mftTrackId, trackTime, timeRes);
    return static_cast<int64_t>(mftStandaloneTracks.lastIndex());
  }

  int64_t fillMchMidStandaloneTrackTable(int64_t mchCollisionId,
                                         int64_t mchMcCollisionId,
                                         int collisionAssignmentStatus,
                                         int64_t mchTrackId,
                                         float trackTime,
                                         float timeRes)
  {
    mchMidStandaloneTracks(mchMcCollisionId, mchCollisionId, collisionAssignmentStatus, mchTrackId, trackTime, timeRes);
    return static_cast<int64_t>(mchMidStandaloneTracks.lastIndex());
  }

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    fGrpMag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());

    if (fGrpMag != nullptr) {
      o2::base::Propagator::initFieldFromGRP(fGrpMag);
      VarManager::SetMagneticField(fGrpMag->getNominalL3Field());
      VarManager::SetupMuonMagField();
      mBzAtMftCenter = fGrpMag->getNominalL3Field();
    }
    mRunNumber = bc.runNumber();
  }

  template <class TBCS>
  void prepareRun(const TBCS& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    const auto bc = bcs.begin();
    initCCDB(bc);
    VarManager::SetMatchingPlane(fzMatching.value);
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(geoPath);
    }

    if (fSamplingFraction < 1.0) {
      std::random_device rd;
      mGenerator = std::mt19937(rd());
    }
    auto hAcceptedEvents = std::get<std::shared_ptr<TH1>>(registry.add("acceptedEvents", "Accepted events", {HistType::kTH1F, {{2, 0, 2.f, ""}}}));
    hAcceptedEvents->GetXaxis()->SetBinLabel(1, "total");
    hAcceptedEvents->GetXaxis()->SetBinLabel(2, "accepted");

    if (fFillQaHistograms) {
      registry.add("eventSelection/vtxZAfterSelection", "Collision vertex z after event selection;V_{z} (cm)", {HistType::kTH1F, {{400, -20, 20, ""}}}, true);
      registry.add("eventSelection/nSelectedMchMidMuonsPerCollision", "Selected MCH-MID muons per collision;N_{#mu}^{MCH-MID};collisions", {HistType::kTH1F, {{51, -0.5, 50.5, ""}}}, true);
      registry.add("allMchMidMftHelix/nMchMidMuons", "MCH-MID muons with collision", {HistType::kTH1F, {{100, 0, 100, ""}}}, true);
      registry.add("allMchMidMftHelix/nMftTracksPerMchMidMuon", "MFT tracks in the same collision per MCH-MID muon", {HistType::kTH1F, {{500, 0, 500, ""}}}, true);
      registry.add("allMchMidMftHelix/zMchMidAtMatching", "MCH-MID z after helix propagation to matching plane;z (cm)", {HistType::kTH1F, {{1000, -200, 200, ""}}}, true);
      registry.add("allMchMidMftHelix/zMftAtMatching", "MFT z after helix propagation to matching plane;z (cm)", {HistType::kTH1F, {{1000, -200, 200, ""}}}, true);
      registry.add("allMchMidMftHelix/debugPMinusPPrime", "MCH p before propagation - MCH p after propagation to MFT z;p-p' (GeV/c)", {HistType::kTH1F, {{2000, -20, 20, ""}}}, true);
      registry.add("mixedMchMidMftHelix/nMchMidMuons", "MCH-MID muons with mixed MFT collisions", {HistType::kTH1F, {{100, 0, 100, ""}}}, true);
      registry.add("mixedMchMidMftHelix/nMixedMftTracksPerMchMidMuon", "Mixed-event MFT tracks per MCH-MID muon", {HistType::kTH1F, {{1000, 0, 1000, ""}}}, true);
      registry.add("mixedMchMidMftHelix/nMixedCollisionsPerMchMidMuon", "Mixed collisions per MCH-MID muon", {HistType::kTH1F, {{100, 0, 100, ""}}}, true);
      registry.add("mixedMchMidMftHelix/debugPMinusPPrime", "mixed MCH p before propagation - MCH p after propagation to MFT z;p-p' (GeV/c)", {HistType::kTH1F, {{2000, -20, 20, ""}}}, true);
      registry.add("bestTimeCollision/nCompatibleCollisions", "Time-compatible collisions per MCH-MID muon;N compatible collisions", {HistType::kTH1F, {{50, 0, 50, ""}}}, true);
      registry.add("bestTimeCollision/bestTimeResidual", "Best collision normalized time residual;#Delta t / #sigma_{t}", {HistType::kTH1F, {{400, -20, 20, ""}}}, true);
      registry.add("bestTimeCollision/secondBestTimeResidual", "Second-best collision normalized time residual;#Delta t / #sigma_{t}", {HistType::kTH1F, {{400, -20, 20, ""}}}, true);
      registry.add("bestTimeCollision/nMftTracksInBestCollision", "MFT tracks in selected best-time collision;N MFT tracks", {HistType::kTH1F, {{500, 0, 500, ""}}}, true);
      registry.add("bestTimeCollision/hasMftInBestCollision", "Whether selected best-time collision has MFT tracks;has MFT", {HistType::kTH1F, {{2, 0, 2, ""}}}, true);
      registry.add("bestTimeCollision/hasMftInOtherCompatibleCollision", "Whether non-best compatible collisions have MFT tracks;has MFT outside best", {HistType::kTH1F, {{2, 0, 2, ""}}}, true);
      registry.add("bestTimeCollision/mftNCompatibleCollisions", "Time-compatible collisions per MFT track;N compatible collisions", {HistType::kTH1F, {{50, 0, 50, ""}}}, true);
      registry.add("bestTimeCollision/mftBestTimeResidual", "MFT best collision normalized time residual;#Delta t / #sigma_{t}", {HistType::kTH1F, {{400, -20, 20, ""}}}, true);
      registry.add("bestTimeCollision/mftSecondBestTimeResidual", "MFT second-best collision normalized time residual;#Delta t / #sigma_{t}", {HistType::kTH1F, {{400, -20, 20, ""}}}, true);
    }

  }

  template <bool isMc, bool withMFTCov, class TCOLLS, class TBCS, class TMUONS, class TMFTS, class TCOVS>
  void fillMchMidMftHelixDeltas(TCOLLS const& collisions,
                                 TBCS const& /*bcs*/,
                                 TMUONS const& muonTracks,
                                 TMFTS const& mftTracks,
                                 TCOVS const& mftCovs)
  {
    std::unordered_set<int64_t> selectedCollisionIds;
    for (const auto& collision : collisions) {
      if (!isSelectedCollision(collision)) {
        continue;
      }
      if (fFillQaHistograms) {
        registry.fill(HIST("eventSelection/vtxZAfterSelection"), collision.posZ());
      }
      const auto collisionId = static_cast<int64_t>(collision.globalIndex());
      selectedCollisionIds.insert(collisionId);
    }

    const bool fillSamePairs = fFillAllMchMidMftHelixDeltas;
    const bool fillMixedPairs = fFillMixedMchMidMftHelixDeltas && fEventMixingDepth.value > 0;
    const bool needMftCache = fillSamePairs || fillMixedPairs;
    const bool needMftTrackCounts = fUseBestTimeCollision.value && fFillQaHistograms.value;

    std::unordered_map<int64_t, std::vector<MftTrackForPairing>> mftsPerCollision;
    std::unordered_map<int64_t, int> mftTrackCountsPerCollision;
    std::unordered_map<int64_t, std::unordered_set<int64_t>> mftMcParticleIdsPerCollision;
    std::unordered_map<int64_t, std::unordered_set<int64_t>> mftMcParticleIdsPerMcCollision;
    std::unordered_map<int64_t, std::vector<PropagatedMchMidTrackAtMatching>> propagatedMchMidMuonsPerCollision;

    if (needMftCache || needMftTrackCounts) {
      // Build the MFT bookkeeping per collision. The expensive propagated cache is only needed for pair output.
      for (const auto& mftTrack : mftTracks) {
        if constexpr (isMc) {
          if (mftTrack.has_mcParticle()) {
            const auto& mftMcParticle = mftTrack.mcParticle();
            if (mftMcParticle.mcCollisionId() >= 0) {
              mftMcParticleIdsPerMcCollision[mftMcParticle.mcCollisionId()].insert(mftMcParticle.globalIndex());
            }
          }
        }

        const auto bestMftCollisionChoice = choosePairCollision<isMc, TBCS>(mftTrack, collisions, selectedCollisionIds);
        if (!bestMftCollisionChoice) {
          continue;
        }
        const auto& bestMftCollision = *bestMftCollisionChoice;
        const int64_t selectedMftCollisionId = bestMftCollision.selectedCollisionId;
        if (selectedCollisionIds.find(selectedMftCollisionId) == selectedCollisionIds.end()) {
          continue;
        }

        if (fUseBestTimeCollision.value && fFillQaHistograms) {
          registry.fill(HIST("bestTimeCollision/mftNCompatibleCollisions"), bestMftCollision.nCompatibleCollisions);
          if (std::isfinite(bestMftCollision.bestTimeResidual)) {
            registry.fill(HIST("bestTimeCollision/mftBestTimeResidual"), bestMftCollision.bestTimeResidual);
          }
          if (std::isfinite(bestMftCollision.secondBestTimeResidual)) {
            registry.fill(HIST("bestTimeCollision/mftSecondBestTimeResidual"), bestMftCollision.secondBestTimeResidual);
          }
        }

        if (needMftTrackCounts) {
          ++mftTrackCountsPerCollision[selectedMftCollisionId];
        }

        if (!needMftCache) {
          continue;
        }

        o2::track::TrackParCovFwd mftTrackPar;
        if constexpr (withMFTCov) {
          auto mftCovIndexIt = mftCovIndexes.find(mftTrack.globalIndex());
          if (mftCovIndexIt == mftCovIndexes.end()) {
            continue;
          }
          auto const& mftTrackCov = mftCovs.rawIteratorAt(mftCovIndexIt->second);
          mftTrackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(mftTrack, 0.f, mftTrackCov);
        } else {
          mftTrackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(mftTrack, 0.f);
        }

        bool mftHasMcParticle = false;
        int64_t mftMcParticleId = -1;
        McTrackOriginInfo mftOrigin;
        if constexpr (isMc) {
          if (mftTrack.has_mcParticle()) {
            mftOrigin = getMcTrackOriginInfo(mftTrack);
            const auto& mftMcParticle = mftTrack.mcParticle();
            mftHasMcParticle = true;
            mftMcParticleId = mftMcParticle.globalIndex();
            mftMcParticleIdsPerCollision[selectedMftCollisionId].insert(mftMcParticleId);
          }
        }

        const int64_t mftTrackId = static_cast<int64_t>(mftTrack.globalIndex());
        const int64_t mftCollisionId = static_cast<int64_t>(mftTrack.collisionId());
        int mftCollisionAssignmentStatus = -1;
        if constexpr (isMc) {
          if (mftCollisionId >= 0 && mftOrigin.mcCollisionId >= 0) {
            const auto& assignedCollision = collisions.rawIteratorAt(mftCollisionId);
            if (assignedCollision.has_mcCollision()) {
              mftCollisionAssignmentStatus = assignedCollision.mcCollisionId() == mftOrigin.mcCollisionId ? 1 : 0;
            }
          }
        }
        const int64_t mftStandaloneTrackId = fillMftStandaloneTrackTable(mftCollisionId,
                                                                           mftOrigin.mcCollisionId,
                                                                           mftCollisionAssignmentStatus,
                                                                           mftTrackId,
                                                                           static_cast<float>(mftTrack.trackTime()),
                                                                           static_cast<float>(mftTrack.trackTimeRes()));
        mftsPerCollision[selectedMftCollisionId].push_back({mftTrackPar,
                                                            makePropagatedTrackValues(mftTrackPar),
                                                            mftTrackId,
                                                            mftStandaloneTrackId,
                                                            absoluteTrackTime<TBCS>(mftTrack, collisions),
                                                            static_cast<float>(mftTrack.trackTimeRes()),
                                                            static_cast<float>(mftTrack.chi2()),
                                                            static_cast<int>(mftTrack.nClusters()),
                                                            mftHasMcParticle,
                                                            mftMcParticleId,
                                                            mftOrigin});
      }
    }

    int nMchMidMuons = 0;
    int nMixedMchMidMuons = 0;
    std::unordered_map<int64_t, int> nSelectedMchMidMuonsPerCollision;
    {
      // Count every MCH-MID track that passes the matching-plane selection before emitting pair rows.
      // Same-event pairs are written while looping over tracks, so the propagated cache alone is incomplete then.
      for (const auto& muonTrack : muonTracks) {
        if (static_cast<int>(muonTrack.trackType()) != 3 ||
            (!isMc && !hasRequiredMchCompatibleCollisionMultiplicity(muonTrack))) {
          continue;
        }
        const auto bestTimeCollisionChoice = choosePairCollision<isMc, TBCS>(muonTrack, collisions, selectedCollisionIds);
        if (!bestTimeCollisionChoice) {
          continue;
        }
        const auto& bestTimeCollision = *bestTimeCollisionChoice;
        const int64_t selectedMuonCollisionId = bestTimeCollision.selectedCollisionId;
        if (selectedCollisionIds.find(selectedMuonCollisionId) == selectedCollisionIds.end()) {
          continue;
        }
        auto muonRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(muonTrack, 0.f, muonTrack);
        auto muonAtMatching = propagateMchToZWithEnergyLoss(muonRaw, fzMatching.value);
        if (!muonAtMatching) {
          continue;
        }
        const auto muonValues = makePropagatedTrackValues(*muonAtMatching);
        if (!isSelectedMchMidAtMatching(muonValues, muonTrack.rAtAbsorberEnd())) {
          continue;
        }
        ++nSelectedMchMidMuonsPerCollision[selectedMuonCollisionId];
      }
    }
    if (fFillQaHistograms) {
      for (const auto selectedCollisionId : selectedCollisionIds) {
        const auto nSelectedMuonIt = nSelectedMchMidMuonsPerCollision.find(selectedCollisionId);
        const int nSelectedMuons = nSelectedMuonIt != nSelectedMchMidMuonsPerCollision.end() ? nSelectedMuonIt->second : 0;
        registry.fill(HIST("eventSelection/nSelectedMchMidMuonsPerCollision"), nSelectedMuons);
      }
    }
    for (const auto& muonTrack : muonTracks) {
      if (static_cast<int>(muonTrack.trackType()) != 3) {
        continue;
      }
      if constexpr (!isMc) {
        if (!hasRequiredMchCompatibleCollisionMultiplicity(muonTrack)) {
          continue;
        }
      }
      const auto bestTimeCollisionChoice = choosePairCollision<isMc, TBCS>(muonTrack, collisions, selectedCollisionIds);
      if (!bestTimeCollisionChoice) {
        continue;
      }
      const auto& bestTimeCollision = *bestTimeCollisionChoice;
      const int64_t selectedMuonCollisionId = bestTimeCollision.selectedCollisionId;
      if (selectedCollisionIds.find(selectedMuonCollisionId) == selectedCollisionIds.end()) {
        continue;
      }

      if (fUseBestTimeCollision.value && fFillQaHistograms) {
        registry.fill(HIST("bestTimeCollision/nCompatibleCollisions"), bestTimeCollision.nCompatibleCollisions);
        if (std::isfinite(bestTimeCollision.bestTimeResidual)) {
          registry.fill(HIST("bestTimeCollision/bestTimeResidual"), bestTimeCollision.bestTimeResidual);
        }
        if (std::isfinite(bestTimeCollision.secondBestTimeResidual)) {
          registry.fill(HIST("bestTimeCollision/secondBestTimeResidual"), bestTimeCollision.secondBestTimeResidual);
        }
        const auto nMftInBestIt = mftTrackCountsPerCollision.find(selectedMuonCollisionId);
        const int nMftInBestCollision = nMftInBestIt != mftTrackCountsPerCollision.end() ? nMftInBestIt->second : 0;
        registry.fill(HIST("bestTimeCollision/nMftTracksInBestCollision"), nMftInBestCollision);
        registry.fill(HIST("bestTimeCollision/hasMftInBestCollision"), nMftInBestCollision > 0 ? 1.f : 0.f);

        bool hasMftInOtherCompatibleCollision = false;
        for (const auto compatibleCollisionId : muonTrack.compatibleCollIds()) {
          const int64_t collisionId = static_cast<int64_t>(compatibleCollisionId);
          if (collisionId == selectedMuonCollisionId) {
            continue;
          }
          const auto nMftIt = mftTrackCountsPerCollision.find(collisionId);
          if (nMftIt != mftTrackCountsPerCollision.end() && nMftIt->second > 0) {
            hasMftInOtherCompatibleCollision = true;
            break;
          }
        }
        registry.fill(HIST("bestTimeCollision/hasMftInOtherCompatibleCollision"), hasMftInOtherCompatibleCollision ? 1.f : 0.f);
      }

      auto muonRaw = o2::aod::fwdtrackutils::getTrackParCovFwdShift(muonTrack, 0.f, muonTrack);
      // The MCH-MID track itself is propagated once to the matching plane and cached for same- and mixed-event pairing.
      auto muonAtMatching = propagateMchToZWithEnergyLoss(muonRaw, fzMatching.value);
      if (!muonAtMatching) {
        continue;
      }
      const auto muonValues = makePropagatedTrackValues(*muonAtMatching);
      if (fFillQaHistograms) {
        registry.fill(HIST("allMchMidMftHelix/zMchMidAtMatching"), muonValues.z);
      }
      if (!isSelectedMchMidAtMatching(muonValues, muonTrack.rAtAbsorberEnd())) {
        continue;
      }
      const auto nSelectedMchMidMuonIt = nSelectedMchMidMuonsPerCollision.find(selectedMuonCollisionId);
      const int nSelectedMchMidMuonsInCollision = nSelectedMchMidMuonIt != nSelectedMchMidMuonsPerCollision.end() ? nSelectedMchMidMuonIt->second : 0;
      McTrackOriginInfo mchOrigin;
      if constexpr (isMc) {
        mchOrigin = getMcTrackOriginInfo(muonTrack);
      }
      const int64_t mchCollisionId = static_cast<int64_t>(muonTrack.collisionId());
      int mchCollisionAssignmentStatus = -1;
      if constexpr (isMc) {
        if (mchCollisionId >= 0 && mchOrigin.mcCollisionId >= 0) {
          const auto& assignedCollision = collisions.rawIteratorAt(mchCollisionId);
          if (assignedCollision.has_mcCollision()) {
            mchCollisionAssignmentStatus = assignedCollision.mcCollisionId() == mchOrigin.mcCollisionId ? 1 : 0;
          }
        }
      }
      const int64_t mchMidStandaloneTrackId = fillMchMidStandaloneTrackTable(mchCollisionId,
                                                                               mchOrigin.mcCollisionId,
                                                                               mchCollisionAssignmentStatus,
                                                                               static_cast<int64_t>(muonTrack.globalIndex()),
                                                                               static_cast<float>(muonTrack.trackTime()),
                                                                               static_cast<float>(muonTrack.trackTimeRes()));
      int64_t mchMidVertexTrackId = -1;
      const auto& collision = collisions.rawIteratorAt(selectedMuonCollisionId);
      auto muonAtVertex = propagateMchToCollisionVertex(muonRaw, collision);
      if (muonAtVertex) {
        mchMidVertexTrackId = fillMchMidVertexTable(selectedMuonCollisionId,
                                                   static_cast<int64_t>(muonTrack.globalIndex()),
                                                   makePropagatedTrackValues(*muonAtVertex));
      }
      const double timeMCH = absoluteTrackTime<TBCS>(muonTrack, collisions);
      const float timeResMCH = static_cast<float>(muonTrack.trackTimeRes());
      propagatedMchMidMuonsPerCollision[selectedMuonCollisionId].push_back({muonRaw,
                                                                            muonValues,
                                                                            static_cast<int64_t>(muonTrack.globalIndex()),
                                                                            mchMidVertexTrackId,
                                                                            mchMidStandaloneTrackId,
                                                                            timeMCH,
                                                                            timeResMCH,
                                                                            static_cast<float>(muonTrack.chi2())});

      if (fillSamePairs) {
        auto mftTracksInCollision = mftsPerCollision.find(selectedMuonCollisionId);
        if (mftTracksInCollision != mftsPerCollision.end()) {
          ++nMchMidMuons;
          int nSelectedMftTracks = 0;
          for (const auto& mftTrackForPair : mftTracksInCollision->second) {
            auto pairPropagation = propagatePair(muonRaw, mftTrackForPair, false);
            if (!pairPropagation) {
              continue;
            }
            ++nSelectedMftTracks;

            McPairInfo mcPairInfo;
            if constexpr (isMc) {
              mcPairInfo = makeMcPairInfo(mchOrigin,
                                          mftTrackForPair,
                                          selectedMuonCollisionId,
                                          mftMcParticleIdsPerCollision,
                                          mftMcParticleIdsPerMcCollision);
            }
            const auto& muonCollision = collisions.rawIteratorAt(selectedMuonCollisionId);
            auto mftAtPVWithMchMomentum = propagateMftToZWithMchMomentum(mftTrackForPair.track,
                                                                         pairPropagation->mchMomentumAtMftZ,
                                                                         static_cast<float>(muonCollision.posZ()));
            const float phiMftPVWithMchMomentum = mftAtPVWithMchMomentum ? mftAtPVWithMchMomentum->phi : std::numeric_limits<float>::quiet_NaN();
            const float etaMftPVWithMchMomentum = mftAtPVWithMchMomentum ? mftAtPVWithMchMomentum->eta : std::numeric_limits<float>::quiet_NaN();
            fillMftMchPairTable(false,
                                mcPairInfo,
                                selectedMuonCollisionId,
                                nSelectedMchMidMuonsInCollision,
                                selectedMuonCollisionId,
                                static_cast<int64_t>(muonTrack.globalIndex()),
                                mftTrackForPair.globalIndex,
                                mchMidVertexTrackId,
                                mftTrackForPair.mftStandaloneTrackId,
                                mchMidStandaloneTrackId,
                                timeMCH,
                                timeResMCH,
                                mftTrackForPair.time,
                                mftTrackForPair.timeRes,
                                muonValues,
                                pairPropagation->mchValuesAtMftLastDisk,
                                static_cast<float>(muonTrack.chi2()),
                                pairPropagation->mftValuesAtMatching,
                                mftTrackForPair.valuesAtMftLastDisk,
                                mftTrackForPair.chi2,
                                mftTrackForPair.nClusters,
                                phiMftPVWithMchMomentum,
                                etaMftPVWithMchMomentum,
                                pairPropagation->mchMomentumAtMftZ,
                                pairPropagation->mchMomentumLoss);
          }
          if (fFillQaHistograms) {
            registry.fill(HIST("allMchMidMftHelix/nMftTracksPerMchMidMuon"), static_cast<float>(nSelectedMftTracks));
          }
        }
      }
    }

    if (fillMixedPairs) {
      using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ>;
      MixingBinning mixingBinning{{fEventMixingVtxBins}, true};
      std::unordered_map<int64_t, std::vector<int>> nMixedCollisionsPerMuon;
      std::unordered_map<int64_t, std::vector<int>> nMixedMftTracksPerMuon;

      // Mix MCH-MID tracks from one selected collision with MFT tracks from another selected collision in the same z-vertex bin.
      auto fillMixedPair = [&](int64_t muonCollisionId, int64_t mftCollisionId) {
        auto muonsInCollision = propagatedMchMidMuonsPerCollision.find(muonCollisionId);
        if (muonsInCollision == propagatedMchMidMuonsPerCollision.end()) {
          return;
        }
        auto mixedMftTracks = mftsPerCollision.find(mftCollisionId);
        if (mixedMftTracks == mftsPerCollision.end()) {
          return;
        }

        auto& mixedCollisionCounts = nMixedCollisionsPerMuon[muonCollisionId];
        auto& mixedMftTrackCounts = nMixedMftTracksPerMuon[muonCollisionId];
        if (mixedCollisionCounts.empty()) {
          mixedCollisionCounts.resize(muonsInCollision->second.size(), 0);
          mixedMftTrackCounts.resize(muonsInCollision->second.size(), 0);
        }

        for (std::size_t iMuon = 0; iMuon < muonsInCollision->second.size(); ++iMuon) {
          ++mixedCollisionCounts[iMuon];
          const auto& muonTrackAtMatching = muonsInCollision->second[iMuon];
          int nSelectedMixedMftTracks = 0;
          for (const auto& mftTrackForPair : mixedMftTracks->second) {
            auto pairPropagation = propagatePair(muonTrackAtMatching.rawTrack, mftTrackForPair, true);
            if (!pairPropagation) {
              continue;
            }
            ++nSelectedMixedMftTracks;
            const auto& muonCollision = collisions.rawIteratorAt(muonCollisionId);
            auto mftAtPVWithMchMomentum = propagateMftToZWithMchMomentum(mftTrackForPair.track,
                                                                         pairPropagation->mchMomentumAtMftZ,
                                                                         static_cast<float>(muonCollision.posZ()));
            const float phiMftPVWithMchMomentum = mftAtPVWithMchMomentum ? mftAtPVWithMchMomentum->phi : std::numeric_limits<float>::quiet_NaN();
            const float etaMftPVWithMchMomentum = mftAtPVWithMchMomentum ? mftAtPVWithMchMomentum->eta : std::numeric_limits<float>::quiet_NaN();
            const auto nSelectedMchMidMuonIt = nSelectedMchMidMuonsPerCollision.find(muonCollisionId);
            const int nSelectedMchMidMuonsInCollision = nSelectedMchMidMuonIt != nSelectedMchMidMuonsPerCollision.end() ? nSelectedMchMidMuonIt->second : 0;
            fillMftMchPairTable(true,
                                McPairInfo{},
                                muonCollisionId,
                                nSelectedMchMidMuonsInCollision,
                                mftCollisionId,
                                muonTrackAtMatching.globalIndex,
                                mftTrackForPair.globalIndex,
                                muonTrackAtMatching.mchMidVertexTrackId,
                                mftTrackForPair.mftStandaloneTrackId,
                                muonTrackAtMatching.mchMidStandaloneTrackId,
                                muonTrackAtMatching.time,
                                muonTrackAtMatching.timeRes,
                                mftTrackForPair.time,
                                mftTrackForPair.timeRes,
                                muonTrackAtMatching.values,
                                pairPropagation->mchValuesAtMftLastDisk,
                                muonTrackAtMatching.chi2,
                                pairPropagation->mftValuesAtMatching,
                                mftTrackForPair.valuesAtMftLastDisk,
                                mftTrackForPair.chi2,
                                mftTrackForPair.nClusters,
                                phiMftPVWithMchMomentum,
                                etaMftPVWithMchMomentum,
                                pairPropagation->mchMomentumAtMftZ,
                                pairPropagation->mchMomentumLoss);
          }
          mixedMftTrackCounts[iMuon] += nSelectedMixedMftTracks;
        }
      };

      for (const auto& [collision1, collision2] : soa::selfCombinations(mixingBinning, fEventMixingDepth.value, -1, collisions, collisions)) {
        if (!isSelectedCollision(collision1) || !isSelectedCollision(collision2)) {
          continue;
        }

        const auto collisionId1 = static_cast<int64_t>(collision1.globalIndex());
        const auto collisionId2 = static_cast<int64_t>(collision2.globalIndex());
        if (collisionId1 == collisionId2) {
          continue;
        }

        if (fMinMixingBcGap.value > 0) {
          const auto bc1 = collision1.template bc_as<TBCS>();
          const auto bc2 = collision2.template bc_as<TBCS>();
          const auto bcGap = std::abs(static_cast<int64_t>(bc1.globalBC()) - static_cast<int64_t>(bc2.globalBC()));
          if (bcGap < fMinMixingBcGap.value) {
            continue;
          }
        }

        fillMixedPair(collisionId1, collisionId2);
        fillMixedPair(collisionId2, collisionId1);
      }

      for (const auto& [collisionId, mixedCollisionCounts] : nMixedCollisionsPerMuon) {
        const auto mftTrackCountsIt = nMixedMftTracksPerMuon.find(collisionId);
        if (mftTrackCountsIt == nMixedMftTracksPerMuon.end()) {
          continue;
        }
        for (std::size_t iMuon = 0; iMuon < mixedCollisionCounts.size(); ++iMuon) {
          if (mixedCollisionCounts[iMuon] <= 0) {
            continue;
          }
          ++nMixedMchMidMuons;
          if (fFillQaHistograms) {
            registry.fill(HIST("mixedMchMidMftHelix/nMixedMftTracksPerMchMidMuon"), mftTrackCountsIt->second[iMuon]);
            registry.fill(HIST("mixedMchMidMftHelix/nMixedCollisionsPerMchMidMuon"), mixedCollisionCounts[iMuon]);
          }
        }
      }
    }

    if (fFillQaHistograms) {
      registry.fill(HIST("allMchMidMftHelix/nMchMidMuons"), nMchMidMuons);
      registry.fill(HIST("mixedMchMidMftHelix/nMchMidMuons"), nMixedMchMidMuons);
    }
  }

  template <bool isMc, bool withMFTCov, class TCOLLS, class TBCS, class TMUONS, class TMFTS, class TCOVS>
  void fillTables(TCOLLS const& collisions,
                  TBCS const& bcs,
                  TMUONS const& muonTracks,
                  TMFTS const& mftTracks,
                  TCOVS const& mftCovs)
  {
    registry.get<TH1>(HIST("acceptedEvents"))->Fill(0);
    if (fSamplingFraction < 1.0 && mDistribution(mGenerator) > fSamplingFraction) {
      return;
    }
    registry.get<TH1>(HIST("acceptedEvents"))->Fill(1);

    if constexpr (withMFTCov) {
      mftCovIndexes.clear();
      for (const auto& mftTrackCov : mftCovs) {
        mftCovIndexes[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
      }
    }

    fillMchMidMftHelixDeltas<isMc, withMFTCov>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  void processMC(MyEventsMC const& collisions,
                 aod::BCsWithTimestamps const& bcs,
                 MyMuonsMC const& muonTracks,
                 MyMFTsMC const& mftTracks,
                 MyMFTCovariances const& mftCovs,
                 aod::McParticles const& /*mcParticles*/)
  {
    prepareRun(bcs);
    fillTables<true, true>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  PROCESS_SWITCH(mftMchPairProducer, processMC, "process_MC", true);

  void processMCNoMFTCov(MyEventsMC const& collisions,
                         aod::BCsWithTimestamps const& bcs,
                         MyMuonsMC const& muonTracks,
                         MyMFTsMC const& mftTracks,
                         aod::McParticles const& /*mcParticles*/)
  {
    prepareRun(bcs);
    fillTables<true, false>(collisions, bcs, muonTracks, mftTracks, nullptr);
  }

  PROCESS_SWITCH(mftMchPairProducer, processMCNoMFTCov, "process_MC_NoMFTCov", false);

  void processRD(MyEvents const& collisions,
                 aod::BCsWithTimestamps const& bcs,
                 MyMuonsWithCov const& muonTracks,
                 MyMFTs const& mftTracks,
                 MyMFTCovariances const& mftCovs)
  {
    prepareRun(bcs);
    fillTables<false, true>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  PROCESS_SWITCH(mftMchPairProducer, processRD, "process_RD", false);

  void processRDNoMFTCov(MyEvents const& collisions,
                         aod::BCsWithTimestamps const& bcs,
                         MyMuonsWithCov const& muonTracks,
                         MyMFTs const& mftTracks)
  {
    prepareRun(bcs);
    fillTables<false, false>(collisions, bcs, muonTracks, mftTracks, nullptr);
  }

  PROCESS_SWITCH(mftMchPairProducer, processRDNoMFTCov, "process_RD_NoMFTCov", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftMchPairProducer>(context)};
};
