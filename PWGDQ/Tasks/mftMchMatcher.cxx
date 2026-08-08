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
/// \file mftMchMatcher.cxx
/// \brief MFT-MCH matching tool for data preparation

#include "PWGDQ/Core/VarManager.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tasks/mftMchMatchingCategories.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
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
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::MFTMults>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMFTs = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels, aod::MFTTrkCompColls>;

namespace matchingcategory = o2::analysis::mftmchmatching;

namespace o2::aod
{
namespace fwdmatchcandidates
{
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TanlMCH, tanlMCH, float);
DECLARE_SOA_COLUMN(InvQPtMCH, invQPtMCH, float);
DECLARE_SOA_COLUMN(TimeMCH, timeMCH, float);
DECLARE_SOA_COLUMN(TimeResMCH, timeResMCH, float);
DECLARE_SOA_COLUMN(Chi2MCH, chi2MCH, float);
DECLARE_SOA_COLUMN(Rabs, rabs, float);
DECLARE_SOA_COLUMN(PDCA, pDCA, float);
DECLARE_SOA_COLUMN(McMaskMCH, mcMaskMCH, int);

DECLARE_SOA_COLUMN(CXXMCH, cXXMCH, float);
DECLARE_SOA_COLUMN(CYYMCH, cYYMCH, float);
DECLARE_SOA_COLUMN(CPhiPhiMCH, cPhiPhiMCH, float);
DECLARE_SOA_COLUMN(CTglTglMCH, cTglTglMCH, float);
DECLARE_SOA_COLUMN(C1Pt1PtMCH, c1Pt1PtMCH, float);
DECLARE_SOA_COLUMN(CXYMCH, cXYMCH, float);
DECLARE_SOA_COLUMN(CPhiXMCH, cPhiXMCH, float);
DECLARE_SOA_COLUMN(CPhiYMCH, cPhiYMCH, float);
DECLARE_SOA_COLUMN(CTglXMCH, cTglXMCH, float);
DECLARE_SOA_COLUMN(CTglYMCH, cTglYMCH, float);
DECLARE_SOA_COLUMN(CTglPhiMCH, cTglPhiMCH, float);
DECLARE_SOA_COLUMN(C1PtXMCH, c1PtXMCH, float);
DECLARE_SOA_COLUMN(C1PtYMCH, c1PtYMCH, float);
DECLARE_SOA_COLUMN(C1PtPhiMCH, c1PtPhiMCH, float);
DECLARE_SOA_COLUMN(C1PtTglMCH, c1PtTglMCH, float);

DECLARE_SOA_COLUMN(XMFT, xMFT, float);
DECLARE_SOA_COLUMN(YMFT, yMFT, float);
DECLARE_SOA_COLUMN(PhiMFT, phiMFT, float);
DECLARE_SOA_COLUMN(TanlMFT, tanlMFT, float);
DECLARE_SOA_COLUMN(InvQPtMFT, invQPtMFT, float);
DECLARE_SOA_COLUMN(TimeMFT, timeMFT, float);
DECLARE_SOA_COLUMN(TimeResMFT, timeResMFT, float);
DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);
DECLARE_SOA_COLUMN(McMaskMFT, mcMaskMFT, int);
DECLARE_SOA_COLUMN(MftClusterSizesAndTrackFlags, mftClusterSizesAndTrackFlags, uint64_t);
DECLARE_SOA_COLUMN(TrackTypeMFT, trackTypeMFT, int);

DECLARE_SOA_COLUMN(CXXMFT, cXXMFT, float);
DECLARE_SOA_COLUMN(CYYMFT, cYYMFT, float);
DECLARE_SOA_COLUMN(CPhiPhiMFT, cPhiPhiMFT, float);
DECLARE_SOA_COLUMN(CTglTglMFT, cTglTglMFT, float);
DECLARE_SOA_COLUMN(C1Pt1PtMFT, c1Pt1PtMFT, float);
DECLARE_SOA_COLUMN(CXYMFT, cXYMFT, float);
DECLARE_SOA_COLUMN(CPhiYMFT, cPhiYMFT, float);
DECLARE_SOA_COLUMN(CPhiXMFT, cPhiXMFT, float);
DECLARE_SOA_COLUMN(CTglXMFT, cTglXMFT, float);
DECLARE_SOA_COLUMN(CTglYMFT, cTglYMFT, float);
DECLARE_SOA_COLUMN(CTglPhiMFT, cTglPhiMFT, float);
DECLARE_SOA_COLUMN(C1PtXMFT, c1PtXMFT, float);
DECLARE_SOA_COLUMN(C1PtYMFT, c1PtYMFT, float);
DECLARE_SOA_COLUMN(C1PtPhiMFT, c1PtPhiMFT, float);
DECLARE_SOA_COLUMN(C1PtTglMFT, c1PtTglMFT, float);

DECLARE_SOA_COLUMN(Chi2Glob, chi2Glob, float);
DECLARE_SOA_COLUMN(Chi2Match, chi2Match, float);
DECLARE_SOA_COLUMN(IsAmbig, isAmbig, bool);
DECLARE_SOA_COLUMN(MFTMult, mftMult, int);
DECLARE_SOA_COLUMN(MatchAttempts, matchAttempts, int);
DECLARE_SOA_COLUMN(DCAX, dcaX, float);
DECLARE_SOA_COLUMN(DCAY, dcaY, float);
DECLARE_SOA_COLUMN(McMaskGlob, mcMaskGlob, int);
DECLARE_SOA_COLUMN(MatchLabel, matchLabel, int);
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);
} // namespace fwdmatchcandidates

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
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
DECLARE_SOA_COLUMN(NSelectedMchMidMuonsInCollision, nSelectedMchMidMuonsInCollision, int);
DECLARE_SOA_COLUMN(MftCollisionId, mftCollisionId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
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

DECLARE_SOA_TABLE(FwdMatchMLCandidates, "AOD", "FWDMLCAND",
                  fwdmatchcandidates::XMCH,
                  fwdmatchcandidates::YMCH,
                  fwdmatchcandidates::PhiMCH,
                  fwdmatchcandidates::TanlMCH,
                  fwdmatchcandidates::InvQPtMCH,
                  fwdmatchcandidates::TimeMCH,
                  fwdmatchcandidates::TimeResMCH,
                  fwdmatchcandidates::Chi2MCH,
                  fwdmatchcandidates::PDCA,
                  fwdmatchcandidates::Rabs,
                  fwdmatchcandidates::CXXMCH,
                  fwdmatchcandidates::CYYMCH,
                  fwdmatchcandidates::CPhiPhiMCH,
                  fwdmatchcandidates::CTglTglMCH,
                  fwdmatchcandidates::C1Pt1PtMCH,
                  fwdmatchcandidates::CXYMCH,
                  fwdmatchcandidates::CPhiYMCH,
                  fwdmatchcandidates::CPhiXMCH,
                  fwdmatchcandidates::CTglXMCH,
                  fwdmatchcandidates::CTglYMCH,
                  fwdmatchcandidates::CTglPhiMCH,
                  fwdmatchcandidates::C1PtXMCH,
                  fwdmatchcandidates::C1PtYMCH,
                  fwdmatchcandidates::C1PtPhiMCH,
                  fwdmatchcandidates::C1PtTglMCH,
                  fwdmatchcandidates::XMFT,
                  fwdmatchcandidates::YMFT,
                  fwdmatchcandidates::PhiMFT,
                  fwdmatchcandidates::TanlMFT,
                  fwdmatchcandidates::InvQPtMFT,
                  fwdmatchcandidates::TimeMFT,
                  fwdmatchcandidates::TimeResMFT,
                  fwdmatchcandidates::Chi2MFT,
                  fwdmatchcandidates::MftClusterSizesAndTrackFlags,
                  fwdmatchcandidates::TrackTypeMFT,
                  fwdmatchcandidates::CXXMFT,
                  fwdmatchcandidates::CYYMFT,
                  fwdmatchcandidates::CPhiPhiMFT,
                  fwdmatchcandidates::CTglTglMFT,
                  fwdmatchcandidates::C1Pt1PtMFT,
                  fwdmatchcandidates::CXYMFT,
                  fwdmatchcandidates::CPhiYMFT,
                  fwdmatchcandidates::CPhiXMFT,
                  fwdmatchcandidates::CTglXMFT,
                  fwdmatchcandidates::CTglYMFT,
                  fwdmatchcandidates::CTglPhiMFT,
                  fwdmatchcandidates::C1PtXMFT,
                  fwdmatchcandidates::C1PtYMFT,
                  fwdmatchcandidates::C1PtPhiMFT,
                  fwdmatchcandidates::C1PtTglMFT,
                  fwdmatchcandidates::Chi2Glob,
                  fwdmatchcandidates::Chi2Match,
                  fwdmatchcandidates::DCAX,
                  fwdmatchcandidates::DCAY,
                  fwdmatchcandidates::IsAmbig,
                  fwdmatchcandidates::MFTMult,
                  fwdmatchcandidates::MatchAttempts,
                  fwdmatchcandidates::McMaskMCH,
                  fwdmatchcandidates::McMaskMFT,
                  fwdmatchcandidates::McMaskGlob,
                  fwdmatchcandidates::MatchLabel,
                  fwdmatchcandidates::IsSignal);

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
                  mftmchmatchpairs::MchCollisionId,
                  mftmchmatchpairs::NSelectedMchMidMuonsInCollision,
                  mftmchmatchpairs::MftCollisionId,
                  mftmchmatchpairs::MchTrackId,
                  mftmchmatchpairs::MftTrackId,
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

DECLARE_SOA_TABLE(MchMidVertexTracks, "AOD", "MCHMIDVTX",
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
>>>>>>> 562aa64bd (daily commit 20260808-1)
} // namespace o2::aod

struct mftMchMatcher {
  Produces<o2::aod::FwdMatchMLCandidates> fwdMatchMLCandidates;

  Produces<o2::aod::MftMchMatchPairs> mftMchMatchPairs;
  Produces<o2::aod::MchMidVertexTracks> mchMidVertexTracks;
  ////   Variables for selecting muon tracks
  Configurable<float> fPMchLow{"cfgPMchLow", 0.0f, ""};
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMFTLow{"cfgEtaMFTlow", -3.6f, ""};
  Configurable<float> fEtaMFTUp{"cfgEtaMFTup", -2.5f, ""};
  Configurable<float> fTrackChi2MFTUp{"cfgTrackChi2MFTUp", 10.f, ""};
  Configurable<float> fPtMFTLow{"cfgPtMFTLow", 0.1f, ""};

  ////   Variables for matching configuration
  Configurable<int> fMaxCandidates{"cfgMaxCandidates", 0, ""};

  Configurable<bool> fKeepBestMatch{"cfgKeepBestMatch", false, "Keep only the best match global muons in the skimming"};
  Configurable<float> fzMatching{"cfgzMatching", -77.5f, "Plane for MFT-MCH matching"};
  Configurable<bool> fFillFwdMatchMLCandidates{"cfgFillFwdMatchMLCandidates", true, "Fill the existing FWDMLCAND table for global MFT-MCH-MID candidates"};
  Configurable<bool> fFillMftMchPairTable{"cfgFillMftMchPairTable", false, "Fill pair-level MCH-MID/MFT matching-plane parameters for same-event and mixed-event studies"};
  Configurable<bool> fFillMchMidVertexTable{"cfgFillMchMidVertexTable", false, "Fill MCH-MID tracks propagated to the collision vertex for dimuon studies"};
  Configurable<bool> fFillAllMchMidMftHelixDeltas{"cfgFillAllMchMidMftHelixDeltas", true, "Process all same-collision MCH-MID/MFT pairs at the matching plane for MFTMCHPAIR output"};
  Configurable<bool> fFillMixedMchMidMftHelixDeltas{"cfgFillMixedMchMidMftHelixDeltas", true, "Process mixed-event MCH-MID/MFT pairs at the matching plane for MFTMCHPAIR output"};
  Configurable<bool> fFillQaHistograms{"cfgFillQaHistograms", false, "Fill QA/debug histograms not needed by the pair-table based downstream analysis"};
  Configurable<bool> fUseBestTimeCollision{"cfgUseBestTimeCollision", false, "Use the time-compatible collision with the smallest normalized time residual for MCH-MID tracks"};
  Configurable<bool> fRequireSingleMchCompatibleCollision{"cfgRequireSingleMchCompatibleCollision", false, "Analyze only MCH-MID tracks with exactly one time-compatible collision"};
  Configurable<int> fEventMixingDepth{"cfgEventMixingDepth", 5, "Number of collisions per vertex-z bin used for event mixing"};
  Configurable<int64_t> fMinMixingBcGap{"cfgMinMixingBcGap", 0, "Minimum absolute global-BC gap between collisions used for MCH-MID/MFT event mixing"};
  Configurable<int> fEventMixingGap{"cfgEventMixingGap", 5, "Deprecated for O2 selfCombinations-based event mixing"};
  ConfigurableAxis fEventMixingVtxBins{"cfgEventMixingVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Vertex-z bins for event mixing"};

  Configurable<float> fSamplingFraction{"cfgSamplingFraction", 1.f, "Fraction of randomly selected events to be processed"};
  Configurable<float> fSamplingFractionTrueLeadingMatches{"cfgSamplingFractionTrueLeadingMatches", 1.f, "Fraction of randomly selected leading true matches to be processed"};
  Configurable<float> fSamplingFractionWrongLeadingMatches{"cfgSamplingFractionWrongLeadingMatches", 1.f, "Fraction of randomly selected leading wrong matches to be processed"};
  Configurable<float> fSamplingFractionDecayLeadingMatches{"cfgSamplingFractionDecayLeadingMatches", 1.f, "Fraction of randomly selected leading decay matches to be processed"};
  Configurable<float> fSamplingFractionFakeLeadingMatches{"cfgSamplingFractionFakeLeadingMatches", 1.f, "Fraction of randomly selected leading fake matches to be processed"};
  Configurable<float> fSamplingFractionTrueNonLeadingMatches{"cfgSamplingFractionTrueNonLeadingMatches", 1.f, "Fraction of randomly selected non-leading true matches to be processed"};
  Configurable<float> fSamplingFractionWrongNonLeadingMatches{"cfgSamplingFractionWrongNonLeadingMatches", 1.f, "Fraction of randomly selected non-leading wrong matches to be processed"};
  Configurable<float> fSamplingFractionDecayNonLeadingMatches{"cfgSamplingFractionDecayNonLeadingMatches", 1.f, "Fraction of randomly selected non-leading decay matches to be processed"};
  Configurable<float> fSamplingFractionFakeNonLeadingMatches{"cfgSamplingFractionFakeNonLeadingMatches", 1.f, "Fraction of randomly selected non-leading fake matches to be processed"};
  Configurable<int> fSamplingBcOddness{"cfgSamplingBcOddness", -1, "Select only events with even (0) or odd (1) global BCs"};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  enum MuonMatchType {
    kMatchTypeTrueLeading = 0,
    kMatchTypeWrongLeading = 1,
    kMatchTypeDecayLeading = 2,
    kMatchTypeFakeLeading = 3,
    kMatchTypeTrueNonLeading = 4,
    kMatchTypeWrongNonLeading = 5,
    kMatchTypeDecayNonLeading = 6,
    kMatchTypeFakeNonLeading = 7,
    kMatchTypeUndefined = 8
  };

  o2::globaltracking::MatchGlobalFwd mExtrap;
  double mBzAtMftCenter{0};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field
  Service<o2::ccdb::BasicCCDBManager> ccdbManager{};
  o2::ccdb::CcdbApi fCCDBApi;

  o2::parameters::GRPMagField* fGrpMag = nullptr;

  std::uniform_real_distribution<double> mDistribution{0.0, 1.0};
  std::mt19937 mGenerator;

  o2::globaltracking::MatchGlobalFwd mMatching;

  std::unordered_map<int64_t, int32_t> mftCovIndexes;

  std::unordered_set<uint32_t> fBestMatch;

  HistogramRegistry registry{"registry", {}};
  
  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack trackToGlobalFwd(const T& track, const C& cov)
  {
    double chi2 = track.chi2();
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> v1{cov.cXX(), cov.cXY(), cov.cYY(), cov.cPhiX(), cov.cPhiY(),
                           cov.cPhiPhi(), cov.cTglX(), cov.cTglY(), cov.cTglPhi(), cov.cTglTgl(),
                           cov.c1PtX(), cov.c1PtY(), cov.c1PtPhi(), cov.c1PtTgl(), cov.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack fwdtrack;
    fwdtrack.setParameters(trackparCov.getParameters());
    fwdtrack.setCovariances(trackparCov.getCovariances());
    return fwdtrack;
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
  };

  struct BestTimeCollisionChoice {
    int64_t selectedCollisionId = -1;
    int64_t originalCollisionId = -1;
    int nCompatibleCollisions = 0;
    float bestTimeResidual = std::numeric_limits<float>::quiet_NaN();
    float secondBestTimeResidual = std::numeric_limits<float>::quiet_NaN();
  };

  using CovarianceValues = std::array<float, 15>;

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

  template <class TCOVARIANCE>
  static CovarianceValues flattenCovariance(const TCOVARIANCE& covariance)
  {
    return {static_cast<float>(covariance(0, 0)),
            static_cast<float>(covariance(1, 1)),
            static_cast<float>(covariance(2, 2)),
            static_cast<float>(covariance(3, 3)),
            static_cast<float>(covariance(4, 4)),
            static_cast<float>(covariance(1, 0)),
            static_cast<float>(covariance(2, 1)),
            static_cast<float>(covariance(2, 0)),
            static_cast<float>(covariance(3, 0)),
            static_cast<float>(covariance(3, 1)),
            static_cast<float>(covariance(3, 2)),
            static_cast<float>(covariance(4, 0)),
            static_cast<float>(covariance(4, 1)),
            static_cast<float>(covariance(4, 2)),
            static_cast<float>(covariance(4, 3))};
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
    if (!fFillMftMchPairTable) {
      return;
    }

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
      mchCollisionId,
      nSelectedMchMidMuonsInCollision,
      mftCollisionId,
      mchTrackId,
      mftTrackId,
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

  void fillMchMidVertexTable(int64_t mchCollisionId,
                             int64_t mchTrackId,
                             const PropagatedTrackAtMatching& vertexValues)
  {
    if (!fFillMchMidVertexTable) {
      return;
    }

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
>>>>>>> 562aa64bd (daily commit 20260808-1)
  }

  template <class T, class C>
  bool passesPDcaCut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static constexpr double sigmaPDCA23 = 80.;
    static constexpr double sigmaPDCA310 = 54.;
    static constexpr double relPRes = 0.0004;
    static constexpr double slopeRes = 0.0005;

    double thetaAbs = TMath::ATan(mchTrack.rAtAbsorberEnd() / 505.) * TMath::RadToDeg();

    // propagate muon track to vertex
    auto trackConv = VarManager::FwdToTrackPar(mchTrack, mchTrack);
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(trackConv.getParameters());
    track.setZ(trackConv.getZ());
    track.setCovariances(trackConv.getCovariances());
    auto mchTrackAtVertex = mMatching.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertex(mchTrackAtVertex, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.p();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = 535. * slopeRes * p;
    double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);

    return (pDCA <= nSigmaPDCA * sigmaPDCAWithRes);
  }

  template <class T, class C>
  bool isGoodMuon(const T& mchTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    // chi2 cut
    if (mchTrack.chi2() > chi2Cut) {
      return false;
    }

    // momentum cut
    if (mchTrack.p() < pCut) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (mchTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mchTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    // RAbs cut
    double rAbs = mchTrack.rAtAbsorberEnd();
    if ((rAbs < rAbsCut[0] || rAbs > rAbsCut[1])) {
      return false;
    }

    // pDCA cut
    if (!passesPDcaCut(mchTrack, collision, nSigmaPdcaCut)) {
      return false;
    }

    return true;
  }

  template <class T>
  bool isGoodMft(const T& mftTrack,
                 double chi2Cut,
                 double pTCut,
                 std::array<double, 2> etaCut)
  {
    // chi2 cut
    if (mftTrack.chi2() > chi2Cut) {
      return false;
    }

    // transverse momentum cut
    if (mftTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mftTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    return true;
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

  static void setMatchTypeLabels(TH1& histogram)
  {
    static constexpr std::array<const char*, 9> labels = {
      "true (leading)",
      "wrong (leading)",
      "decay (leading)",
      "fake (leading)",
      "true (non leading)",
      "wrong (non leading)",
      "decay (non leading)",
      "fake (non leading)",
      "undefined"};
    for (std::size_t i = 0; i < labels.size(); ++i) {
      histogram.GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, labels[i]);
    }
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    fCCDBApi.init(ccdburl);
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

    AxisSpec matchTypeAxis = {static_cast<int>(kMatchTypeUndefined) + 1, 0, static_cast<double>(kMatchTypeUndefined) + 1, ""};
    auto hMatchType = std::get<std::shared_ptr<TH1>>(registry.add("matchType", "Match type", {HistType::kTH1F, {matchTypeAxis}}));
    hMatchType->GetXaxis()->SetBinLabel(1, "true (leading)");
    hMatchType->GetXaxis()->SetBinLabel(2, "wrong (leading)");
    hMatchType->GetXaxis()->SetBinLabel(3, "decay (leading)");
    hMatchType->GetXaxis()->SetBinLabel(4, "fake (leading)");
    hMatchType->GetXaxis()->SetBinLabel(5, "true (non leading)");
    hMatchType->GetXaxis()->SetBinLabel(6, "wrong (non leading)");
    hMatchType->GetXaxis()->SetBinLabel(7, "decay (non leading)");
    hMatchType->GetXaxis()->SetBinLabel(8, "fake (non leading)");
    hMatchType->GetXaxis()->SetBinLabel(9, "undefined");

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

    if (fFillFwdMatchMLCandidates && fFillQaHistograms) {
      AxisSpec matchTypeAxis = {static_cast<int>(kMatchTypeUndefined) + 1, 0, static_cast<double>(kMatchTypeUndefined) + 1, ""};
      auto hMatchType = std::get<std::shared_ptr<TH1>>(registry.add("matchType", "Match type", {HistType::kTH1F, {matchTypeAxis}}));
      setMatchTypeLabels(*hMatchType);
>>>>>>> 562aa64bd (daily commit 20260808-1)

      auto hMatchTypeAccepted = std::get<std::shared_ptr<TH1>>(registry.add("matchTypeAccepted", "Match type (accepted)", {HistType::kTH1F, {matchTypeAxis}}));
      setMatchTypeLabels(*hMatchTypeAccepted);
    }
  }

  template <typename TMuons>
  void fillBestMuonMatches(TMuons const& muons)
  {
    fBestMatch.clear();
    std::unordered_map<int, std::pair<float, int>> mCandidates;
    for (const auto& muon : muons) {
      if (static_cast<int>(muon.trackType()) < 2) {
        const auto muonID = muon.matchMCHTrackId();
        const auto chi2 = muon.chi2MatchMCHMFT();
        auto [candidate, inserted] = mCandidates.try_emplace(muonID, chi2, muon.globalIndex());
        if (!inserted && chi2 < candidate->second.first) {
          candidate->second = {chi2, muon.globalIndex()};
        }
      }
    }
    for (const auto& pairCand : mCandidates) {
      fBestMatch[pairCand.second.second] = true;
=======
    for (const auto& candidateByMuon : mCandidates) {
      fBestMatch.emplace(candidateByMuon.second.second);
>>>>>>> 562aa64bd (daily commit 20260808-1)
    }
  }

  template <class TMUON, class TMFT>
  void fillMatchablePairs(TMUON const& muonTracks,
                          TMFT const& mftTracks,
                          std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    // outer loop on muon tracks
    for (const auto& muonTrack : muonTracks) {
      // only consider MCH standalone or MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) <= 2) {
        continue;
      }

      // only consider tracks associated to some collision
      if (!muonTrack.has_collision()) {
        continue;
      }

      // skip tracks that do not have an associated MC particle
      if (!muonTrack.has_mcParticle()) {
        continue;
      }
      // get the index associated to the MC particle
      auto muonMcParticle = muonTrack.mcParticle();

      int64_t muonMcTrackIndex = muonMcParticle.globalIndex();

      // inner loop on MFT tracks
      for (const auto& mftTrack : mftTracks) {
        // skip tracks that do not have an associated MC particle
        if (!mftTrack.has_mcParticle()) {
          continue;
        }
        // get the index associated to the MC particle
        auto mftMcParticle = mftTrack.mcParticle();
        int64_t mftMcTrackIndex = mftMcParticle.globalIndex();

        if (muonMcTrackIndex == mftMcTrackIndex) {
          matchablePairs.emplace_back(static_cast<int64_t>(muonTrack.globalIndex()),
                                      static_cast<int64_t>(mftTrack.globalIndex()));
        }
      }
    }
  }

  static bool isPairedMuon(int64_t muonTrackId, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    for (const auto& [id1, id2] : matchablePairs) {
      if (muonTrackId == id1) {
        return true;
      }
    }
    return false;
    return std::any_of(matchablePairs.begin(), matchablePairs.end(), [muonTrackId](const auto& pair) {
      return muonTrackId == pair.first;
    });
  template <class EVT, class BC, class TMUON, class TMFTS>
  int getMftMchMatchAttempts(EVT const& collisions,
                             BC const& bcs,
                             TMUON const& mchTrack,
                             TMFTS const& mftTracks)
  {
    if (!mchTrack.has_collision()) {
      return 0;
    }
    const auto& collMch = collisions.rawIteratorAt(mchTrack.collisionId());
    const auto& bcMch = bcs.rawIteratorAt(collMch.bcId());

    int attempts{0};
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_collision()) {
        continue;
      }

      const auto& collMft = collisions.rawIteratorAt(mftTrack.collisionId());
      const auto& bcMft = bcs.rawIteratorAt(collMft.bcId());

      int64_t deltaBc = static_cast<int64_t>(bcMft.globalBC()) - static_cast<int64_t>(bcMch.globalBC());
      double deltaBcNS = o2::constants::lhc::LHCBunchSpacingNS * deltaBc;
      double deltaTrackTime = mftTrack.trackTime() - mchTrack.trackTime() + deltaBcNS;
      double trackTimeResTot = mftTrack.trackTimeRes() + mchTrack.trackTimeRes();

      if (std::fabs(deltaTrackTime) > trackTimeResTot) {
        continue;
      }
      attempts += 1;
    }

    return attempts;
=======

  float getMatchTypeSamplingFraction(MuonMatchType matchType) const
  {
    switch (matchType) {
      case kMatchTypeTrueLeading:
        return fSamplingFractionTrueLeadingMatches.value;
      case kMatchTypeWrongLeading:
        return fSamplingFractionWrongLeadingMatches.value;
      case kMatchTypeDecayLeading:
        return fSamplingFractionDecayLeadingMatches.value;
      case kMatchTypeFakeLeading:
        return fSamplingFractionFakeLeadingMatches.value;
      case kMatchTypeTrueNonLeading:
        return fSamplingFractionTrueNonLeadingMatches.value;
      case kMatchTypeWrongNonLeading:
        return fSamplingFractionWrongNonLeadingMatches.value;
      case kMatchTypeDecayNonLeading:
        return fSamplingFractionDecayNonLeadingMatches.value;
      case kMatchTypeFakeNonLeading:
        return fSamplingFractionFakeNonLeadingMatches.value;
      default:
        return 1.f;
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

    const bool fillSamePairs = fFillAllMchMidMftHelixDeltas && fFillMftMchPairTable;
    const bool fillMixedPairs = fFillMixedMchMidMftHelixDeltas && fFillMftMchPairTable && fEventMixingDepth.value > 0;
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

        if (!mftTrack.has_collision()) {
          continue;
        }
        const auto bestMftCollision = chooseBestTimeCollision<TBCS>(mftTrack, collisions, selectedCollisionIds);
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

        mftsPerCollision[selectedMftCollisionId].push_back({mftTrackPar,
                                                            makePropagatedTrackValues(mftTrackPar),
                                                            static_cast<int64_t>(mftTrack.globalIndex()),
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
    if (fFillMftMchPairTable || fFillQaHistograms) {
      // Count every MCH-MID track that passes the matching-plane selection before emitting pair rows.
      // Same-event pairs are written while looping over tracks, so the propagated cache alone is incomplete then.
      for (const auto& muonTrack : muonTracks) {
        if (static_cast<int>(muonTrack.trackType()) != 3 || !muonTrack.has_collision() ||
            !hasRequiredMchCompatibleCollisionMultiplicity(muonTrack)) {
          continue;
        }
        const auto bestTimeCollision = chooseBestTimeCollision<TBCS>(muonTrack, collisions, selectedCollisionIds);
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
      if (!muonTrack.has_collision()) {
        continue;
      }
      if (!hasRequiredMchCompatibleCollisionMultiplicity(muonTrack)) {
        continue;
      }
      const auto bestTimeCollision = chooseBestTimeCollision<TBCS>(muonTrack, collisions, selectedCollisionIds);
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
      propagatedMchMidMuonsPerCollision[selectedMuonCollisionId].push_back({muonRaw, muonValues, static_cast<int64_t>(muonTrack.globalIndex()), static_cast<float>(muonTrack.chi2())});

      if (fFillMchMidVertexTable) {
        const auto& collision = collisions.rawIteratorAt(selectedMuonCollisionId);
        auto muonAtVertex = propagateMchToCollisionVertex(muonRaw, collision);
        if (muonAtVertex) {
          fillMchMidVertexTable(selectedMuonCollisionId,
                                static_cast<int64_t>(muonTrack.globalIndex()),
                                makePropagatedTrackValues(*muonAtVertex));
        }
      }

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

  template <bool isMc, class TCOLLS, class TBCS, class TMUONS, class TMFTS, class TCOVS>
  void fillTable(TCOLLS const& collisions,
                 TBCS const& bcs,
                 TMUONS const& muonTracks,
                 TMFTS const& mftTracks,
                 TCOVS const& mftCovs)
  {
    std::unordered_map<int64_t, int> matchAttemptsMap;

    registry.get<TH1>(HIST("acceptedEvents"))->Fill(0);
    // reject a randomly selected fraction of events
    if (fSamplingFraction < 1.0) {
      double rnd = mDistribution(mGenerator);
      if (rnd > fSamplingFraction) {
        return;
      }
    }
    registry.get<TH1>(HIST("acceptedEvents"))->Fill(1);

    if (fFillFwdMatchMLCandidates && (fKeepBestMatch.value || isMc)) {
      fillBestMuonMatches(muonTracks);
    }

    std::vector<std::pair<int64_t, int64_t>> matchablePairs;
    if constexpr (isMc) {
      if (fFillFwdMatchMLCandidates) {
        fillMatchablePairs(muonTracks, mftTracks, matchablePairs);
      }
    }

    mftCovIndexes.clear();
    for (const auto& mftTrackCov : mftCovs) {
      mftCovIndexes[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    for (const auto& muon : muonTracks) {
    if (((fFillAllMchMidMftHelixDeltas || fFillMixedMchMidMftHelixDeltas) && fFillMftMchPairTable) || fFillMchMidVertexTable) {
      fillMchMidMftHelixDeltas<isMc, withMFTCov>(collisions, bcs, muonTracks, mftTracks, mftCovs);
    }

    if (!fFillFwdMatchMLCandidates) {
      return;
    }

    for (auto muon : muonTracks) {
      // only consider global MFT-MCH-MID matches
      if (static_cast<int>(muon.trackType()) != 0) {
        continue;
      }

      if (!muon.has_collision()) {
        continue;
      }

      bool isBestMatch = fBestMatch.find(muon.globalIndex()) != fBestMatch.end();

      if (fKeepBestMatch && !isBestMatch) {
        continue;
      }

      const auto& collision = collisions.rawIteratorAt(muon.collisionId());
      auto bc_coll = collision.template bc_as<TBCS>();

      auto muontrack = muon.template matchMCHTrack_as<TMUONS>();
      if (!hasRequiredMchCompatibleCollisionMultiplicity(muontrack)) {
        continue;
      }
      auto mfttrack = muon.template matchMFTTrack_as<TMFTS>();
      auto const& mfttrackcov = mftCovs.rawIteratorAt(mftCovIndexes[mfttrack.globalIndex()]);

      const auto muonTime = muontrack.trackTime() + bc_coll.globalBC() * o2::constants::lhc::LHCBunchSpacingNS;
      const auto mftTime = mfttrack.trackTime() + bc_coll.globalBC() * o2::constants::lhc::LHCBunchSpacingNS;

      o2::track::TrackParCovFwd mftprop = VarManager::FwdToTrackPar(mfttrack, mfttrackcov);
      o2::dataformats::GlobalFwdTrack muonprop = trackToGlobalFwd(muontrack, muontrack);
      if (fzMatching.value < 0.) {
        mftprop = VarManager::PropagateFwd(mfttrack, mfttrackcov, fzMatching.value);
        muonprop = VarManager::PropagateMuon(muontrack, collision, VarManager::kToMatching);
      }
      auto muonpropCov = muonprop.getCovariances();
      auto mftpropCov = mftprop.getCovariances();
      const auto muonCovValues = flattenCovariance(muonprop.getCovariances());

      float mftX = mfttrack.x();
      float mftY = mfttrack.y();
      float mftPhi = mfttrack.phi();
      float mftTgl = mfttrack.tgl();
      float mftInvQPt = mfttrack.signed1Pt();
      CovarianceValues mftCovValues;
      mftCovValues.fill(-999.f);

      if constexpr (withMFTCov) {
        auto mftCovIndexIt = mftCovIndexes.find(mfttrack.globalIndex());
        if (mftCovIndexIt == mftCovIndexes.end()) {
          continue;
        }
        auto const& mfttrackcov = mftCovs.rawIteratorAt(mftCovIndexIt->second);
        o2::track::TrackParCovFwd mftprop = VarManager::FwdToTrackPar(mfttrack, mfttrackcov);
        if (fzMatching.value < 0.) {
          mftprop = VarManager::PropagateFwd(mfttrack, mfttrackcov, fzMatching.value);
        }
        mftX = mftprop.getX();
        mftY = mftprop.getY();
        mftPhi = mftprop.getPhi();
        mftTgl = mftprop.getTgl();
        mftInvQPt = mftprop.getInvQPt();
        mftCovValues = flattenCovariance(mftprop.getCovariances());
      }

      if (!isGoodMuon(muontrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMFTLow, fEtaMFTUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp)) {
        continue;
      }

      // at this level we consider all the matching candidates, regardless of the MFT tracks quality
      // MFT track quality cuts should be applied only after having selected the best candidate
      // if (!isGoodMft(mfttrack, fTrackChi2MFTUp, fPtMFTLow, {fEtaMFTLow, fEtaMFTUp})) {
      //  continue;
      // }

      bool IsAmbig = (muon.compatibleCollIds().size() != 1);
      int MFTMult = collision.mftNtracks();
      int matchAttempts = 0;
      auto matchAttemptsIt = matchAttemptsMap.find(muontrack.globalIndex());
      if (matchAttemptsIt == matchAttemptsMap.end()) {
        matchAttempts = getMftMchMatchAttempts(collisions, bcs, muontrack, mftTracks);
        matchAttemptsMap.insert(std::make_pair(static_cast<int64_t>(muontrack.globalIndex()), matchAttempts));
      } else {
        matchAttempts = matchAttemptsIt->second;
      }
      const bool isAmbiguous = muon.compatibleCollIds().size() != 1;
      constexpr int mftMultiplicity = -1;

      auto matchType = kMatchTypeUndefined;
      if constexpr (isMc) {
        matchType = getMatchType(muon, muonTracks, mftTracks, matchablePairs, isBestMatch);
      }
      const bool isSignal = (matchType == kMatchTypeTrueLeading) || (matchType == kMatchTypeTrueNonLeading);

      int mcMaskMuon = 0;
      int mcMaskMft = 0;
      int mcMaskGlob = 0;
      if constexpr (isMc) {
        mcMaskMuon = muontrack.mcMask();
        mcMaskMft = mfttrack.mcMask();
        mcMaskGlob = muon.mcMask();
      }

      if (fFillQaHistograms) {
        registry.get<TH1>(HIST("matchType"))->Fill(static_cast<int>(matchType));
      }

      // skipp odd/even BCs if requested
      if (fSamplingBcOddness.value >= 0 && (static_cast<int>((bc_coll.globalBC() % 2)) != fSamplingBcOddness.value)) {
        continue;
      }

      const float matchTypeSamplingFraction = getMatchTypeSamplingFraction(matchType);

      if (matchTypeSamplingFraction < 1.0) {
        double rnd = mDistribution(mGenerator);
        if (rnd > matchTypeSamplingFraction) {
          continue;
        }
      }

      if (fFillQaHistograms) {
        registry.get<TH1>(HIST("matchTypeAccepted"))->Fill(static_cast<int>(matchType));
      }

      fwdMatchMLCandidates(
        muonprop.getX(),
        muonprop.getY(),
        muonprop.getPhi(),
        muonprop.getTgl(),
        muonprop.getInvQPt(),
        muonTime,
        muontrack.trackTimeRes(),
        muontrack.chi2(),
        muontrack.pDca(),
        muontrack.rAtAbsorberEnd(),
        muonpropCov(0, 0),
        muonpropCov(1, 1),
        muonpropCov(2, 2),
        muonpropCov(3, 3),
        muonpropCov(4, 4),
        muonpropCov(1, 0),
        muonpropCov(2, 1),
        muonpropCov(2, 0),
        muonpropCov(3, 0),
        muonpropCov(3, 1),
        muonpropCov(3, 2),
        muonpropCov(4, 0),
        muonpropCov(4, 1),
        muonpropCov(4, 2),
        muonpropCov(4, 3),
        mftprop.getX(),
        mftprop.getY(),
        mftprop.getPhi(),
        mftprop.getTgl(),
        mftprop.getInvQPt(),
        muonCovValues[0],
        muonCovValues[1],
        muonCovValues[2],
        muonCovValues[3],
        muonCovValues[4],
        muonCovValues[5],
        muonCovValues[6],
        muonCovValues[7],
        muonCovValues[8],
        muonCovValues[9],
        muonCovValues[10],
        muonCovValues[11],
        muonCovValues[12],
        muonCovValues[13],
        muonCovValues[14],
        mftX,
        mftY,
        mftPhi,
        mftTgl,
        mftInvQPt,
        mftTime,
        mfttrack.trackTimeRes(),
        mfttrack.chi2(),
        mfttrack.mftClusterSizesAndTrackFlags(),
        (mfttrack.isCA() ? 1 : 0),
        mftpropCov(0, 0),
        mftpropCov(1, 1),
        mftpropCov(2, 2),
        mftpropCov(3, 3),
        mftpropCov(4, 4),
        mftpropCov(1, 0),
        mftpropCov(2, 1),
        mftpropCov(2, 0),
        mftpropCov(3, 0),
        mftpropCov(3, 1),
        mftpropCov(3, 2),
        mftpropCov(4, 0),
        mftpropCov(4, 1),
        mftpropCov(4, 2),
        mftpropCov(4, 3),
        muon.chi2(),
        muon.chi2MatchMCHMFT(),
        muon.fwdDcaX(),
        muon.fwdDcaY(),
        IsAmbig,
        MFTMult,
        matchAttempts,
        isAmbiguous,
        mftMultiplicity,
        mcMaskMuon,
        mcMaskMft,
        mcMaskGlob,
        static_cast<int>(matchType),
        isSignal);
    }
  }

  void processMC(MyEvents const& collisions,
                 aod::BCsWithTimestamps const& bcs,
                 MyMuonsMC const& muonTracks,
                 MyMFTsMC const& mftTracks,
                 MyMFTCovariances const& mftCovs,
                 aod::McParticles const& /*mcParticles*/)
  {
    if (bcs.size() > 0) {
      auto bc = bcs.begin();
      initCCDB(bc);
      VarManager::SetMatchingPlane(fzMatching.value);
    }

    fillTable<true>(collisions, bcs, muonTracks, mftTracks, mftCovs);
    prepareRun(bcs);
    fillTable<true, true>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  PROCESS_SWITCH(mftMchMatcher, processMC, "process_MC", true);

  void processMCNoMFTCov(MyEvents const& collisions,
                         aod::BCsWithTimestamps const& bcs,
                         MyMuonsMC const& muonTracks,
                         MyMFTsMC const& mftTracks,
                         aod::McParticles const& /*mcParticles*/)
  {
    prepareRun(bcs);
    fillTable<true, false>(collisions, bcs, muonTracks, mftTracks, nullptr);
  }

  PROCESS_SWITCH(mftMchMatcher, processMCNoMFTCov, "process_MC_NoMFTCov", false);

  void processRD(MyEvents const& collisions,
                 aod::BCsWithTimestamps const& bcs,
                 MyMuonsWithCov const& muonTracks,
                 MyMFTs const& mftTracks,
                 MyMFTCovariances const& mftCovs)
  {
    if (bcs.size() > 0) {
      auto bc = bcs.begin();
      initCCDB(bc);
      VarManager::SetMatchingPlane(fzMatching.value);
    }

    fillTable<false>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  PROCESS_SWITCH(mftMchMatcher, processRD, "process_RD", false);
=======
    prepareRun(bcs);
    fillTable<false, true>(collisions, bcs, muonTracks, mftTracks, mftCovs);
  }

  PROCESS_SWITCH(mftMchMatcher, processRD, "process_RD", false);

  void processRDNoMFTCov(MyEvents const& collisions,
                         aod::BCsWithTimestamps const& bcs,
                         MyMuonsWithCov const& muonTracks,
                         MyMFTs const& mftTracks)
  {
    prepareRun(bcs);
    fillTable<false, false>(collisions, bcs, muonTracks, mftTracks, nullptr);
  }

  PROCESS_SWITCH(mftMchMatcher, processRDNoMFTCov, "process_RD_NoMFTCov", false);
>>>>>>> 562aa64bd (daily commit 20260808-1)
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftMchMatcher>(context)};
};
