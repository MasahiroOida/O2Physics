// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file check_globalTrack.cxx
/// \brief Diagnostic study of multiple type-0 global-muon candidates per MCH-MID track.

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct CheckGlobalTrack {
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {21, -0.5, 20.5}, "candidate multiplicity"};
  ConfigurableAxis axisX{"axisX", {400, -100., 100.}, "#it{x} (cm)"};
  ConfigurableAxis axisY{"axisY", {400, -100., 100.}, "#it{y} (cm)"};
  ConfigurableAxis axisEta{"axisEta", {250, -4.5, -2.0}, "#it{#eta}"};
  ConfigurableAxis axisPhi{"axisPhi", {180, -3.141592653589793, 3.141592653589793}, "#it{#varphi} (rad)"};
  ConfigurableAxis axisPt{"axisPt", {200, 0., 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisMatchingChi2{"axisMatchingChi2", {500, 0., 100.}, "MFT-MCH matching #chi^{2}"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct TrackValues {
    float x = 0.f;
    float y = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    float pt = 0.f;
  };

  struct Candidate {
    int64_t globalTrackId = -1;
    int64_t mftTrackId = -1;
    float matchingChi2 = 0.f;
    TrackValues global;
    TrackValues mft;
  };

  struct MchCandidateGroup {
    TrackValues mch;
    std::vector<Candidate> candidates;
  };

  enum class TrackHistogramGroup {
    GlobalCandidates,
    UniqueMch,
    DistinctMft,
    BestGlobalCandidate,
    BestMft
  };

  void addTrackHistograms(const std::string& prefix, const std::string& objectTitle)
  {
    const std::string titleX = objectTitle + ";#it{x} (cm);tracks";
    const std::string titleY = objectTitle + ";#it{y} (cm);tracks";
    const std::string titleEta = objectTitle + ";#it{#eta};tracks";
    const std::string titlePhi = objectTitle + ";#it{#varphi} (rad);tracks";
    const std::string titlePt = objectTitle + ";#it{p}_{T} (GeV/#it{c});tracks";
    registry.add(prefix + "x", titleX.c_str(), HistType::kTH1F, {axisX});
    registry.add(prefix + "y", titleY.c_str(), HistType::kTH1F, {axisY});
    registry.add(prefix + "eta", titleEta.c_str(), HistType::kTH1F, {axisEta});
    registry.add(prefix + "phi", titlePhi.c_str(), HistType::kTH1F, {axisPhi});
    registry.add(prefix + "pt", titlePt.c_str(), HistType::kTH1F, {axisPt});
  }

  void init(InitContext const&)
  {
    auto candidateStatus = registry.add<TH1>("candidateStatus",
                                              "Type-0 candidate validation;selection step;rows",
                                              HistType::kTH1I, {{3, -0.5, 2.5}});
    candidateStatus->GetXaxis()->SetBinLabel(1, "trackType == 0");
    candidateStatus->GetXaxis()->SetBinLabel(2, "valid MCH and MFT relations");
    candidateStatus->GetXaxis()->SetBinLabel(3, "MCH relation is type 3");

    auto bestMatchStatus = registry.add<TH1>("bestMatchStatus",
                                              "Best-match selection per MCH-MID track;selection step;MCH-MID tracks",
                                              HistType::kTH1I, {{3, -0.5, 2.5}});
    bestMatchStatus->GetXaxis()->SetBinLabel(1, "valid type-0 group");
    bestMatchStatus->GetXaxis()->SetBinLabel(2, ">= 2 distinct MFT tracks");
    bestMatchStatus->GetXaxis()->SetBinLabel(3, "positive finite chi2 candidate");

    registry.add("multiplicity/nType0RowsPerMch",
                 "Valid type-0 rows per MCH-MID track;type-0 rows;MCH-MID tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("multiplicity/nDistinctMftTracksPerMch",
                 "Distinct linked MFT tracks per MCH-MID track;distinct MFT tracks;MCH-MID tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("multiplicity/nRepeatedRowsPerMch",
                 "Repeated type-0 rows with an already used MFT relation per MCH-MID track;rows - distinct MFT tracks;MCH-MID tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("multiplicity/nPositiveChi2CandidatesPerMch",
                 "Candidates eligible for best-match selection per MCH-MID track;positive finite matching-#chi^{2} candidates;MCH-MID tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("multiplicity/nDistinctMftVsNType0Rows",
                 "Distinct MFT relations versus type-0 rows per MCH-MID track;type-0 rows;distinct MFT tracks",
                 HistType::kTH2I, {axisMultiplicity, axisMultiplicity});
    registry.add("duplicateChecks/nType0RowsPerMchMftPair",
                 "Valid type-0 rows per reconstructed MCH-MFT pair;type-0 rows per pair;MCH-MFT pairs",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("duplicateChecks/nSelectedBestPerMch",
                 "Selected best global tracks per reconstructed MCH-MID track in the multi-MFT sample;selected global tracks per MCH-MID track;MCH-MID tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("duplicateChecks/nSelectedBestPerMft",
                 "Selected best global tracks per reconstructed MFT track in the multi-MFT sample;selected global tracks per MFT track;MFT tracks",
                 HistType::kTH1I, {axisMultiplicity});
    registry.add("duplicateChecks/nSelectedBestPerMchMftPair",
                 "Selected best rows per reconstructed MCH-MFT pair in the multi-MFT sample;selected rows per pair;MCH-MFT pairs",
                 HistType::kTH1I, {axisMultiplicity});

    addTrackHistograms("multipleMft/globalCandidates/",
                       "Type-0 global candidates in MCH groups with multiple distinct MFT tracks");
    addTrackHistograms("multipleMft/uniqueMch/",
                       "Unique MCH-MID tracks with multiple distinct MFT candidates");
    addTrackHistograms("multipleMft/distinctMft/",
                       "Distinct MFT tracks linked to the same MCH-MID track");
    addTrackHistograms("multipleMft/bestGlobalCandidate/",
                       "Minimum matching-chi2 type-0 global candidate per MCH-MID track");
    addTrackHistograms("multipleMft/bestMft/",
                       "MFT track of the minimum matching-chi2 candidate per MCH-MID track");
    registry.add("multipleMft/globalCandidates/matchingChi2",
                 "Type-0 candidates in MCH groups with multiple distinct MFT tracks;MFT-MCH matching #chi^{2};global candidates",
                 HistType::kTH1F, {axisMatchingChi2});
    registry.add("multipleMft/bestGlobalCandidate/matchingChi2",
                 "Minimum positive MFT-MCH matching #chi^{2} per MCH-MID track;MFT-MCH matching #chi^{2};selected global candidates",
                 HistType::kTH1F, {axisMatchingChi2});
  }

  static float wrappedPhi(float phi)
  {
    return std::atan2(std::sin(phi), std::cos(phi));
  }

  template <typename TTrack>
  static TrackValues values(const TTrack& track)
  {
    return {track.x(), track.y(), track.eta(), wrappedPhi(track.phi()), track.pt()};
  }

  void fillTrackHistograms(TrackHistogramGroup group, const TrackValues& track)
  {
#define FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS(prefix) \
  registry.fill(HIST(prefix "x"), track.x);         \
  registry.fill(HIST(prefix "y"), track.y);         \
  registry.fill(HIST(prefix "eta"), track.eta);     \
  registry.fill(HIST(prefix "phi"), track.phi);     \
  registry.fill(HIST(prefix "pt"), track.pt)

    switch (group) {
      case TrackHistogramGroup::GlobalCandidates:
        FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS("multipleMft/globalCandidates/");
        break;
      case TrackHistogramGroup::UniqueMch:
        FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS("multipleMft/uniqueMch/");
        break;
      case TrackHistogramGroup::DistinctMft:
        FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS("multipleMft/distinctMft/");
        break;
      case TrackHistogramGroup::BestGlobalCandidate:
        FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS("multipleMft/bestGlobalCandidate/");
        break;
      case TrackHistogramGroup::BestMft:
        FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS("multipleMft/bestMft/");
        break;
    }

#undef FILL_CHECK_GLOBAL_TRACK_HISTOGRAMS
  }

  void process(aod::FwdTracks const& fwdTracks, aod::MFTTracks const& /*mftTracks*/)
  {
    std::unordered_map<int64_t, MchCandidateGroup> groupsByMchTrack;
    std::unordered_map<int64_t, std::unordered_map<int64_t, int>> nType0RowsPerMchMftPair;
    std::unordered_map<int64_t, int> nSelectedBestPerMch;
    std::unordered_map<int64_t, int> nSelectedBestPerMft;
    std::unordered_map<int64_t, std::unordered_map<int64_t, int>> nSelectedBestPerMchMftPair;

    // Keep all type-0 rows here. No event, kinematic, or matching-chi2 cut is applied.
    for (const auto& globalTrack : fwdTracks) {
      if (globalTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        continue;
      }
      registry.fill(HIST("candidateStatus"), 0.);

      if (!globalTrack.has_matchMCHTrack() || !globalTrack.has_matchMFTTrack()) {
        continue;
      }
      registry.fill(HIST("candidateStatus"), 1.);

      const auto mchTrack = globalTrack.template matchMCHTrack_as<aod::FwdTracks>();
      const auto mftTrack = globalTrack.template matchMFTTrack_as<aod::MFTTracks>();
      if (mchTrack.trackType() != aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        continue;
      }
      registry.fill(HIST("candidateStatus"), 2.);

      const int64_t mchTrackId = mchTrack.globalIndex();
      const int64_t mftTrackId = mftTrack.globalIndex();
      auto [groupIterator, inserted] = groupsByMchTrack.try_emplace(mchTrackId);
      if (inserted) {
        groupIterator->second.mch = values(mchTrack);
      }
      ++nType0RowsPerMchMftPair[mchTrackId][mftTrackId];
      groupIterator->second.candidates.push_back({globalTrack.globalIndex(),
                                                  mftTrackId,
                                                  globalTrack.chi2MatchMCHMFT(),
                                                  values(globalTrack),
                                                  values(mftTrack)});
    }

    for (const auto& groupEntry : groupsByMchTrack) {
      const auto& group = groupEntry.second;
      std::unordered_map<int64_t, const Candidate*> distinctMftCandidates;
      const Candidate* bestCandidate = nullptr;
      int nPositiveChi2Candidates = 0;
      for (const auto& candidate : group.candidates) {
        // The MFT parameters are identical for repeated references to the same MFT row.
        distinctMftCandidates.try_emplace(candidate.mftTrackId, &candidate);
        if (!std::isfinite(candidate.matchingChi2) || candidate.matchingChi2 <= 0.f) {
          continue;
        }
        ++nPositiveChi2Candidates;
        // A strict comparison keeps the first candidate when matching chi2 is equal.
        if (bestCandidate == nullptr || candidate.matchingChi2 < bestCandidate->matchingChi2) {
          bestCandidate = &candidate;
        }
      }

      const int nType0Rows = static_cast<int>(group.candidates.size());
      const int nDistinctMftTracks = static_cast<int>(distinctMftCandidates.size());
      registry.fill(HIST("bestMatchStatus"), 0.);
      registry.fill(HIST("multiplicity/nType0RowsPerMch"), nType0Rows);
      registry.fill(HIST("multiplicity/nDistinctMftTracksPerMch"), nDistinctMftTracks);
      registry.fill(HIST("multiplicity/nRepeatedRowsPerMch"), nType0Rows - nDistinctMftTracks);
      registry.fill(HIST("multiplicity/nPositiveChi2CandidatesPerMch"), nPositiveChi2Candidates);
      registry.fill(HIST("multiplicity/nDistinctMftVsNType0Rows"), nType0Rows, nDistinctMftTracks);

      if (nDistinctMftTracks < 2) {
        continue;
      }
      registry.fill(HIST("bestMatchStatus"), 1.);

      fillTrackHistograms(TrackHistogramGroup::UniqueMch, group.mch);
      for (const auto& candidate : group.candidates) {
        fillTrackHistograms(TrackHistogramGroup::GlobalCandidates, candidate.global);
        if (std::isfinite(candidate.matchingChi2)) {
          registry.fill(HIST("multipleMft/globalCandidates/matchingChi2"), candidate.matchingChi2);
        }
      }
      for (const auto& distinctEntry : distinctMftCandidates) {
        fillTrackHistograms(TrackHistogramGroup::DistinctMft, distinctEntry.second->mft);
      }
      if (bestCandidate == nullptr) {
        continue;
      }
      registry.fill(HIST("bestMatchStatus"), 2.);
      fillTrackHistograms(TrackHistogramGroup::BestGlobalCandidate, bestCandidate->global);
      fillTrackHistograms(TrackHistogramGroup::BestMft, bestCandidate->mft);
      registry.fill(HIST("multipleMft/bestGlobalCandidate/matchingChi2"), bestCandidate->matchingChi2);
      ++nSelectedBestPerMch[groupEntry.first];
      ++nSelectedBestPerMft[bestCandidate->mftTrackId];
      ++nSelectedBestPerMchMftPair[groupEntry.first][bestCandidate->mftTrackId];
    }

    for (const auto& mchEntry : nType0RowsPerMchMftPair) {
      for (const auto& pairEntry : mchEntry.second) {
        registry.fill(HIST("duplicateChecks/nType0RowsPerMchMftPair"), pairEntry.second);
      }
    }
    for (const auto& entry : nSelectedBestPerMch) {
      registry.fill(HIST("duplicateChecks/nSelectedBestPerMch"), entry.second);
    }
    for (const auto& entry : nSelectedBestPerMft) {
      registry.fill(HIST("duplicateChecks/nSelectedBestPerMft"), entry.second);
    }
    for (const auto& mchEntry : nSelectedBestPerMchMftPair) {
      for (const auto& pairEntry : mchEntry.second) {
        registry.fill(HIST("duplicateChecks/nSelectedBestPerMchMftPair"), pairEntry.second);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckGlobalTrack>(cfgc)};
}
