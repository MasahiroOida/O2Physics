// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file globalTrackDimuonStudy.cxx
/// \brief Same- and mixed-event dimuon mass spectra from global-track matching candidates.

#include "Common/Tasks/globalTrackMatchingStudyDataModel.h"
#include "Common/Tasks/mftMchMatchingCategories.h"

#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace matchingcategory = o2::analysis::mftmchmatching;

struct GlobalTrackDimuonStudy {
  struct Candidate {
    float signed1Pt = 0.f;
    float pt = 0.f;
    float eta = 0.f;
    float phi = 0.f;
  };

  using CandidatesByCategory = std::array<std::vector<Candidate>, 6>;
  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ>;

  static constexpr int firstStudyCategory = 1;
  static constexpr int lastStudyCategory = 5;

  Configurable<int> cfgMixingDepth{"cfgMixingDepth", 5, "Number of previous collisions mixed in the same z-vertex bin"};
  ConfigurableAxis cfgMixingVtxBins{"cfgMixingVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins in z vertex (cm)"};
  ConfigurableAxis axisMass{"axisMass", {600, 0.f, 6.f}, "#it{m}_{#mu#mu} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPairPt{"axisPairPt", {200, 0.f, 20.f}, "#it{p}_{T}^{#mu#mu} (GeV/#it{c})"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr float massMuon = 0.1056583745f; // GeV/c^2

  static int categoryIndex(int category) { return category + 1; }

  static std::string categoryName(int index)
  {
    using Category = matchingcategory::McMatchStudyCategory;
    switch (static_cast<Category>(index - 1)) {
      case Category::Unknown:
        return "unknown";
      case Category::CorrectPrimaryMuon:
        return "correctPrimaryMuon";
      case Category::CorrectPiKDecayMuon:
        return "correctPiKDecayMuon";
      case Category::CorrectOtherMuon:
        return "correctOtherMuon";
      case Category::FalseWithTrueCounterpart:
        return "falseWithTrueCounterpart";
      case Category::NoTrueMftCounterpart:
        return "noTrueMftCounterpart";
    }
    return "invalid";
  }

  static std::string categoryPairName(int first, int second)
  {
    if (first > second) {
      std::swap(first, second);
    }
    return categoryName(first) + "_" + categoryName(second);
  }

  static bool isValidCandidate(const Candidate& candidate)
  {
    return std::isfinite(candidate.signed1Pt) && std::isfinite(candidate.pt) && std::isfinite(candidate.eta) && std::isfinite(candidate.phi) &&
           candidate.signed1Pt != 0.f && candidate.pt > 0.f;
  }

  static float invariantMass(const Candidate& first, const Candidate& second)
  {
    const float firstPx = first.pt * std::cos(first.phi);
    const float firstPy = first.pt * std::sin(first.phi);
    const float firstPz = first.pt * std::sinh(first.eta);
    const float secondPx = second.pt * std::cos(second.phi);
    const float secondPy = second.pt * std::sin(second.phi);
    const float secondPz = second.pt * std::sinh(second.eta);
    const float firstEnergy = std::sqrt(firstPx * firstPx + firstPy * firstPy + firstPz * firstPz + massMuon * massMuon);
    const float secondEnergy = std::sqrt(secondPx * secondPx + secondPy * secondPy + secondPz * secondPz + massMuon * massMuon);
    return std::sqrt(std::max((firstEnergy + secondEnergy) * (firstEnergy + secondEnergy) -
                              (firstPx + secondPx) * (firstPx + secondPx) -
                              (firstPy + secondPy) * (firstPy + secondPy) -
                              (firstPz + secondPz) * (firstPz + secondPz),
                              0.f));
  }

  void addPairHistograms(const std::string& pairType, const std::string& category)
  {
    const std::string prefix = pairType + "/" + category + "/";
    registry.add(prefix + "ulsMass", "Unlike-sign dimuon invariant mass;#it{m}_{#mu#mu} (GeV/#it{c}^{2});pairs", HistType::kTH1F, {axisMass});
    registry.add(prefix + "lsppMass", "Like-sign ++ dimuon invariant mass;#it{m}_{#mu#mu} (GeV/#it{c}^{2});pairs", HistType::kTH1F, {axisMass});
    registry.add(prefix + "lsmmMass", "Like-sign -- dimuon invariant mass;#it{m}_{#mu#mu} (GeV/#it{c}^{2});pairs", HistType::kTH1F, {axisMass});
    registry.add(prefix + "ulsMassVsPt", "Unlike-sign dimuon invariant mass versus pair transverse momentum;#it{m}_{#mu#mu} (GeV/#it{c}^{2});#it{p}_{T}^{#mu#mu} (GeV/#it{c})", HistType::kTH2F, {axisMass, axisPairPt});
  }

  void init(InitContext const&)
  {
    if (cfgMixingDepth.value < 1) {
      LOGF(fatal, "cfgMixingDepth must be at least one");
    }
    addPairHistograms("same", "all");
    addPairHistograms("mixed", "all");
    for (int first = firstStudyCategory; first <= lastStudyCategory; ++first) {
      for (int second = first; second <= lastStudyCategory; ++second) {
        const auto category = categoryPairName(first, second);
        addPairHistograms("same", category);
        addPairHistograms("mixed", category);
      }
    }
    registry.add("event/nCandidates", "Selected global-track candidates per collision;#it{N}_{tracks};collisions", HistType::kTH1I, {{101, -0.5, 100.5}});
    registry.add("event/nMixedCollisionPairs", "Collision pairs used for event mixing;collision pairs;entries", HistType::kTH1I, {{1, -0.5, 0.5}});
  }

  void fillPair(const Candidate& first, const Candidate& second, bool sameEvent, int firstCategory, int secondCategory, bool fillAll)
  {
    const float mass = invariantMass(first, second);
    const float pairPt = std::hypot(first.pt * std::cos(first.phi) + second.pt * std::cos(second.phi),
                                    first.pt * std::sin(first.phi) + second.pt * std::sin(second.phi));

#define FILL_DIMUON_PAIR(prefix)                                                     \
  if (first.signed1Pt * second.signed1Pt < 0.f) {                                    \
    registry.fill(HIST(prefix "ulsMass"), mass);                                    \
    registry.fill(HIST(prefix "ulsMassVsPt"), mass, pairPt);                        \
  } else if (first.signed1Pt > 0.f) {                                                \
    registry.fill(HIST(prefix "lsppMass"), mass);                                   \
  } else {                                                                           \
    registry.fill(HIST(prefix "lsmmMass"), mass);                                   \
  }

#define FILL_CATEGORY_PAIR(name)                                                     \
    if (sameEvent) {                                                                  \
      FILL_DIMUON_PAIR("same/" name "/")                                            \
    } else {                                                                          \
      FILL_DIMUON_PAIR("mixed/" name "/")                                           \
    }

    if (fillAll) {
      if (sameEvent) {
        FILL_DIMUON_PAIR("same/all/")
      } else {
        FILL_DIMUON_PAIR("mixed/all/")
      }
    }

    if (firstCategory > secondCategory) {
      std::swap(firstCategory, secondCategory);
    }
    switch (firstCategory * 10 + secondCategory) {
      case 11: FILL_CATEGORY_PAIR("correctPrimaryMuon_correctPrimaryMuon"); break;
      case 12: FILL_CATEGORY_PAIR("correctPrimaryMuon_correctPiKDecayMuon"); break;
      case 13: FILL_CATEGORY_PAIR("correctPrimaryMuon_correctOtherMuon"); break;
      case 14: FILL_CATEGORY_PAIR("correctPrimaryMuon_falseWithTrueCounterpart"); break;
      case 15: FILL_CATEGORY_PAIR("correctPrimaryMuon_noTrueMftCounterpart"); break;
      case 22: FILL_CATEGORY_PAIR("correctPiKDecayMuon_correctPiKDecayMuon"); break;
      case 23: FILL_CATEGORY_PAIR("correctPiKDecayMuon_correctOtherMuon"); break;
      case 24: FILL_CATEGORY_PAIR("correctPiKDecayMuon_falseWithTrueCounterpart"); break;
      case 25: FILL_CATEGORY_PAIR("correctPiKDecayMuon_noTrueMftCounterpart"); break;
      case 33: FILL_CATEGORY_PAIR("correctOtherMuon_correctOtherMuon"); break;
      case 34: FILL_CATEGORY_PAIR("correctOtherMuon_falseWithTrueCounterpart"); break;
      case 35: FILL_CATEGORY_PAIR("correctOtherMuon_noTrueMftCounterpart"); break;
      case 44: FILL_CATEGORY_PAIR("falseWithTrueCounterpart_falseWithTrueCounterpart"); break;
      case 45: FILL_CATEGORY_PAIR("falseWithTrueCounterpart_noTrueMftCounterpart"); break;
      case 55: FILL_CATEGORY_PAIR("noTrueMftCounterpart_noTrueMftCounterpart"); break;
      default: break;
    }

#undef FILL_CATEGORY_PAIR
#undef FILL_DIMUON_PAIR
  }

  void fillPairs(const CandidatesByCategory& first, const CandidatesByCategory& second, bool sameEvent)
  {
    if (sameEvent) {
      for (int firstCategory = 0; firstCategory < static_cast<int>(first.size()); ++firstCategory) {
        for (size_t i = 0; i < first[firstCategory].size(); ++i) {
          for (size_t j = i + 1; j < first[firstCategory].size(); ++j) {
            fillPair(first[firstCategory][i], first[firstCategory][j], true, firstCategory, firstCategory, true);
          }
          for (int secondCategory = firstCategory + 1; secondCategory < static_cast<int>(first.size()); ++secondCategory) {
            for (const auto& secondCandidate : first[secondCategory]) {
              fillPair(first[firstCategory][i], secondCandidate, true, firstCategory, secondCategory, true);
            }
          }
        }
      }
    } else {
      for (int firstCategory = 0; firstCategory < static_cast<int>(first.size()); ++firstCategory) {
        for (const auto& firstCandidate : first[firstCategory]) {
          for (int secondCategory = 0; secondCategory < static_cast<int>(second.size()); ++secondCategory) {
            for (const auto& secondCandidate : second[secondCategory]) {
              fillPair(firstCandidate, secondCandidate, false, firstCategory, secondCategory, true);
            }
          }
        }
      }
    }
  }

  void process(aod::Collisions const& collisions, aod::GlobalTrackMatchingStudies const& tracks)
  {
    std::unordered_map<int64_t, CandidatesByCategory> candidatesByCollision;
    candidatesByCollision.reserve(tracks.size());
    for (const auto& track : tracks) {
      const int64_t collisionId = track.collisionId();
      const int index = categoryIndex(track.mcMatchStudyCategory());
      const Candidate candidate{track.signed1PtMCH(), track.ptMCH(), track.etaMCH(), track.phiMCH()};
      if (collisionId < 0 || index < 0 || index >= static_cast<int>(CandidatesByCategory{}.size()) || !isValidCandidate(candidate)) {
        continue;
      }
      candidatesByCollision[collisionId][index].push_back(candidate);
    }

    for (const auto& collision : collisions) {
      const auto candidates = candidatesByCollision.find(collision.globalIndex());
      if (candidates == candidatesByCollision.end()) {
        continue;
      }
      size_t nCandidates = 0;
      for (const auto& category : candidates->second) {
        nCandidates += category.size();
      }
      registry.fill(HIST("event/nCandidates"), nCandidates);
      fillPairs(candidates->second, candidates->second, true);
    }

    MixingBinning mixingBinning{{cfgMixingVtxBins}, true};
    for (const auto& [firstCollision, secondCollision] : soa::selfCombinations(mixingBinning, cfgMixingDepth.value, -1, collisions, collisions)) {
      const int64_t firstCollisionId = firstCollision.globalIndex();
      const int64_t secondCollisionId = secondCollision.globalIndex();
      if (firstCollisionId == secondCollisionId) {
        continue;
      }
      const auto first = candidatesByCollision.find(firstCollisionId);
      const auto second = candidatesByCollision.find(secondCollisionId);
      if (first == candidatesByCollision.end() || second == candidatesByCollision.end()) {
        continue;
      }
      registry.fill(HIST("event/nMixedCollisionPairs"), 0.);
      fillPairs(first->second, second->second, false);
    }
  }
  PROCESS_SWITCH(GlobalTrackDimuonStudy, process, "Build same- and mixed-event dimuon pairs from GlobalTrackMatchingStudies", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalTrackDimuonStudy>(cfgc)};
}
