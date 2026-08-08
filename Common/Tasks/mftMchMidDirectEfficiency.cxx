// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file mftMchMidDirectEfficiency.cxx
/// \brief Direct reconstructed MFT-MCH-MID / MCH-MID track-count study.

#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

using DirectStudyCollisions = soa::Join<aod::Collisions, aod::EvSels>;

namespace o2::aod::mftmchmiddirect
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(ReferencedMchTrackId, referencedMchTrackId, int64_t);
DECLARE_SOA_COLUMN(GlobalTrackId, globalTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(HasMftMchMidTrack, hasMftMchMidTrack, bool);
DECLARE_SOA_COLUMN(NGlobalMatchesForMch, nGlobalMatchesForMch, int);
DECLARE_SOA_COLUMN(ZVtx, zVtx, float);
DECLARE_SOA_COLUMN(RAbs, rAbs, float);
DECLARE_SOA_COLUMN(ZMCH, zMCH, float);
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TglMCH, tglMCH, float);
DECLARE_SOA_COLUMN(PtMCH, ptMCH, float);
DECLARE_SOA_COLUMN(EtaMCH, etaMCH, float);
DECLARE_SOA_COLUMN(PMCH, pMCH, float);
DECLARE_SOA_COLUMN(Chi2MCH, chi2MCH, float);
DECLARE_SOA_COLUMN(Chi2MatchMCHMID, chi2MatchMCHMID, float);
DECLARE_SOA_COLUMN(EtaMFT, etaMFT, float);
} // namespace o2::aod::mftmchmiddirect

namespace o2::aod
{
DECLARE_SOA_TABLE(MftMchMidDirectProbes, "AOD", "MFTMCHDIRECT",
                  mftmchmiddirect::CollisionId,
                  mftmchmiddirect::MchTrackId,
                  mftmchmiddirect::ReferencedMchTrackId,
                  mftmchmiddirect::GlobalTrackId,
                  mftmchmiddirect::MftTrackId,
                  mftmchmiddirect::HasMftMchMidTrack,
                  mftmchmiddirect::NGlobalMatchesForMch,
                  mftmchmiddirect::ZVtx,
                  mftmchmiddirect::RAbs,
                  mftmchmiddirect::ZMCH,
                  mftmchmiddirect::XMCH,
                  mftmchmiddirect::YMCH,
                  mftmchmiddirect::PhiMCH,
                  mftmchmiddirect::TglMCH,
                  mftmchmiddirect::PtMCH,
                  mftmchmiddirect::EtaMCH,
                  mftmchmiddirect::PMCH,
                  mftmchmiddirect::Chi2MCH,
                  mftmchmiddirect::Chi2MatchMCHMID,
                  mftmchmiddirect::EtaMFT);
} // namespace o2::aod

struct MftMchMidDirectEfficiency
{
  Configurable<bool> applySel8{"cfgApplySel8", true, "Require sel8 collisions"};
  Configurable<float> etaMchMin{"cfgEtaMchMin", -3.6f, "Strict lower MCH-MID eta boundary"};
  Configurable<float> etaMchMax{"cfgEtaMchMax", -2.5f, "Strict upper MCH-MID eta boundary"};
  Configurable<float> ptMchMin{"cfgPtMchMin", 0.3f, "Strict lower MCH-MID pT boundary (GeV/c)"};
  Configurable<float> rAbsMin{"cfgRAbsMin", 27.5f, "Strict lower MCH-MID absorber-radius boundary (cm)"};
  Configurable<float> rAbsMax{"cfgRAbsMax", 89.5f, "Strict upper MCH-MID absorber-radius boundary (cm)"};
  Configurable<float> etaMftMin{"cfgEtaMftMin", -3.6f, "Strict lower referenced-MFT eta boundary"};
  Configurable<float> etaMftMax{"cfgEtaMftMax", -2.5f, "Strict upper referenced-MFT eta boundary"};

  Produces<aod::MftMchMidDirectProbes> directProbes;

  struct DirectGlobalMatch {
    int64_t collisionId = -1;
    int64_t globalTrackId = -1;
    int64_t referencedMchTrackId = -1;
    int64_t mftTrackId = -1;
    float etaMft = std::numeric_limits<float>::quiet_NaN();
  };

  template <typename T>
  bool passesMchCuts(const T& track) const
  {
    return std::isfinite(track.pt()) && std::isfinite(track.eta()) &&
           track.pt() > ptMchMin.value &&
           track.eta() > etaMchMin.value && track.eta() < etaMchMax.value &&
           std::isfinite(track.rAtAbsorberEnd()) &&
           track.rAtAbsorberEnd() > rAbsMin.value && track.rAtAbsorberEnd() < rAbsMax.value;
  }

  bool passesMftEta(float eta) const
  {
    return std::isfinite(eta) && eta > etaMftMin.value && eta < etaMftMax.value;
  }

  void process(DirectStudyCollisions const& collisions,
               aod::FwdTracks const& fwdTracks,
               aod::MFTTracks const& mftTracks)
  {
    std::unordered_map<int64_t, float> selectedCollisionZ;
    selectedCollisionZ.reserve(collisions.size());
    for (const auto& collision : collisions) {
      if (!applySel8.value || collision.sel8()) {
        selectedCollisionZ.emplace(collision.globalIndex(), collision.posZ());
      }
    }

    std::unordered_map<int64_t, float> mftEtaById;
    mftEtaById.reserve(mftTracks.size());
    for (const auto& mftTrack : mftTracks) {
      mftEtaById.emplace(mftTrack.globalIndex(), mftTrack.eta());
    }

    std::unordered_map<int64_t, int64_t> mchCollisionById;
    for (const auto& mchTrack : fwdTracks) {
      if (mchTrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack &&
          mchTrack.has_collision()) {
        mchCollisionById.emplace(mchTrack.globalIndex(), mchTrack.collisionId());
      }
    }

    std::unordered_map<int64_t, std::vector<DirectGlobalMatch>> globalMatchesByMch;
    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack ||
          !globalTrack.has_collision() ||
          selectedCollisionZ.find(globalTrack.collisionId()) == selectedCollisionZ.end()) {
        continue;
      }
      const int64_t referencedMchTrackId = globalTrack.matchMCHTrackId();
      const int64_t mftTrackId = globalTrack.matchMFTTrackId();
      const auto mchIt = mchCollisionById.find(referencedMchTrackId);
      const auto mftIt = mftEtaById.find(mftTrackId);
      if (referencedMchTrackId < 0 || mchIt == mchCollisionById.end()) {
        LOGF(fatal, "trackType=0 globalIndex=%lld has invalid matchMCHTrackId=%lld",
             static_cast<long long>(globalTrack.globalIndex()), static_cast<long long>(referencedMchTrackId));
        continue;
      }
      if (mchIt->second != globalTrack.collisionId()) {
        LOGF(fatal, "Collision mismatch for trackType=0 globalIndex=%lld: global collision=%lld, referenced MCH-MID collision=%lld",
             static_cast<long long>(globalTrack.globalIndex()), static_cast<long long>(globalTrack.collisionId()),
             static_cast<long long>(mchIt->second));
        continue;
      }
      if (mftTrackId < 0 || mftIt == mftEtaById.end()) {
        LOGF(fatal, "trackType=0 globalIndex=%lld has invalid matchMFTTrackId=%lld",
             static_cast<long long>(globalTrack.globalIndex()), static_cast<long long>(mftTrackId));
        continue;
      }
      if (!passesMftEta(mftIt->second)) {
        continue;
      }
      globalMatchesByMch[referencedMchTrackId].push_back({globalTrack.collisionId(),
                                                          globalTrack.globalIndex(),
                                                          referencedMchTrackId,
                                                          mftTrackId,
                                                          mftIt->second});
    }
    for (auto& entry : globalMatchesByMch) {
      auto& matches = entry.second;
      std::sort(matches.begin(), matches.end(), [](const auto& first, const auto& second) {
        return first.globalTrackId < second.globalTrackId;
      });
    }

    for (const auto& mchTrack : fwdTracks) {
      if (mchTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack ||
          !mchTrack.has_collision()) {
        continue;
      }
      const auto collisionIt = selectedCollisionZ.find(mchTrack.collisionId());
      if (collisionIt == selectedCollisionZ.end() || !passesMchCuts(mchTrack)) {
        continue;
      }

      const DirectGlobalMatch* selectedMatch = nullptr;
      int nGlobalMatchesForMch = 0;
      const auto matchIt = globalMatchesByMch.find(mchTrack.globalIndex());
      if (matchIt != globalMatchesByMch.end()) {
        for (const auto& match : matchIt->second) {
          if (match.collisionId != mchTrack.collisionId()) {
            continue;
          }
          ++nGlobalMatchesForMch;
          if (!selectedMatch) {
            selectedMatch = &match;
          }
        }
      }

      const bool hasMftMchMidTrack = selectedMatch != nullptr;
      const int64_t invalidId = -1;
      const float invalidFloat = std::numeric_limits<float>::quiet_NaN();
      directProbes(mchTrack.collisionId(),
                   mchTrack.globalIndex(),
                   hasMftMchMidTrack ? selectedMatch->referencedMchTrackId : invalidId,
                   hasMftMchMidTrack ? selectedMatch->globalTrackId : invalidId,
                   hasMftMchMidTrack ? selectedMatch->mftTrackId : invalidId,
                   hasMftMchMidTrack,
                   nGlobalMatchesForMch,
                   collisionIt->second,
                   mchTrack.rAtAbsorberEnd(),
                   mchTrack.z(),
                   mchTrack.x(),
                   mchTrack.y(),
                   mchTrack.phi(),
                   mchTrack.tgl(),
                   mchTrack.pt(),
                   mchTrack.eta(),
                   mchTrack.p(),
                   mchTrack.chi2(),
                   mchTrack.chi2MatchMCHMID(),
                   hasMftMchMidTrack ? selectedMatch->etaMft : invalidFloat);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MftMchMidDirectEfficiency>(cfgc)};
}
