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

/// \file mftTrackPairMass.cxx
/// \brief Same-event MFT standalone track-pair mass study with MC truth labels.

#include "Common/DataModel/CollisionAssociationTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct MFTTrackPairMass {
  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MFTTrackLabeled = MFTTracksLabeled::iterator;

  enum TruthCategory : int {
    kUnknown = 0,
    kPiPi = 1,
    kKK = 2,
    kMuMu = 3,
    kPiK = 4,
    kOther = 5
  };

  enum MassHypothesis : int {
    kPiPiMass = 0,
    kKKMass = 1,
    kMuMuMass = 2
  };

  struct FourVector {
    float px = 0.f;
    float py = 0.f;
    float pz = 0.f;
    float e = 0.f;
  };

  struct PairKinematics {
    float mass = std::numeric_limits<float>::quiet_NaN();
    float pt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float openingAngle = std::numeric_limits<float>::quiet_NaN();
  };

  struct TruthMuonMother {
    bool valid = false;
    int motherPdg = 0;
    int64_t motherId = -1;
  };

  Configurable<float> minTrackPt{"minTrackPt", 0.f, "Minimum MFT track pT"};
  Configurable<float> maxTrackPt{"maxTrackPt", 100.f, "Maximum MFT track pT"};
  Configurable<float> minTrackEta{"minTrackEta", -5.f, "Minimum MFT track eta"};
  Configurable<float> maxTrackEta{"maxTrackEta", -2.f, "Maximum MFT track eta"};
  Configurable<int> minNClusters{"minNClusters", 0, "Minimum number of MFT clusters"};
  Configurable<float> maxTrackChi2{"maxTrackChi2", 1.e9f, "Maximum MFT track chi2"};
  Configurable<int> maxTracksPerCollision{"maxTracksPerCollision", 300, "Maximum selected tracks per collision used for pairing"};
  Configurable<int> nBinsMass{"nBinsMass", 600, "Number of pair mass bins"};
  Configurable<float> minMass{"minMass", 0.f, "Minimum pair mass"};
  Configurable<float> maxMass{"maxMass", 6.f, "Maximum pair mass"};
  Configurable<int> nBinsPairPt{"nBinsPairPt", 200, "Number of pair pT bins"};
  Configurable<float> minPairPt{"minPairPt", 0.f, "Minimum pair pT"};
  Configurable<float> maxPairPt{"maxPairPt", 20.f, "Maximum pair pT"};
  Configurable<int> nBinsPairEta{"nBinsPairEta", 120, "Number of pair eta bins"};
  Configurable<float> minPairEta{"minPairEta", -6.f, "Minimum pair eta"};
  Configurable<float> maxPairEta{"maxPairEta", 0.f, "Maximum pair eta"};
  Configurable<int> nBinsOpeningAngle{"nBinsOpeningAngle", 160, "Number of opening-angle bins"};
  Configurable<float> minOpeningAngle{"minOpeningAngle", 0.f, "Minimum opening angle"};
  Configurable<float> maxOpeningAngle{"maxOpeningAngle", 3.2f, "Maximum opening angle"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMass{nBinsMass.value, minMass.value, maxMass.value, "m (GeV/#it{c}^{2})"};
    AxisSpec axisPairPt{nBinsPairPt.value, minPairPt.value, maxPairPt.value, "p_{T}^{pair} (GeV/#it{c})"};
    AxisSpec axisPairEta{nBinsPairEta.value, minPairEta.value, maxPairEta.value, "#eta^{pair}"};
    AxisSpec axisOpeningAngle{nBinsOpeningAngle.value, minOpeningAngle.value, maxOpeningAngle.value, "#alpha (rad)"};
    AxisSpec axisTruthCategory{6, -0.5, 5.5, "truth category"};

    registry.add("hAssociationMultiplicity", "MFT association multiplicity;N_{coll};tracks", kTH1F, {{10, -0.5, 9.5}});
    registry.add("hSelectedTracksPerCollision", "selected MFT tracks per collision;N_{tracks};collisions", kTH1F, {{500, -0.5, 499.5}});
    registry.add("hSkippedCollisionsByMultiplicity", "collisions skipped by track multiplicity;N_{tracks};collisions", kTH1F, {{500, -0.5, 499.5}});
    registry.add("hTrackPtEta", "selected MFT tracks;p_{T} (GeV/#it{c});#eta", kTH2F, {{200, 0., 20.}, {120, -6., 0.}});

    addPairHistograms("PiPi", axisMass, axisPairPt, axisPairEta, axisOpeningAngle, axisTruthCategory);
    addPairHistograms("KK", axisMass, axisPairPt, axisPairEta, axisOpeningAngle, axisTruthCategory);
    addPairHistograms("MuMu", axisMass, axisPairPt, axisPairEta, axisOpeningAngle, axisTruthCategory);
    registry.add("OS/hPairTruthCategory", "opposite-sign pair truth category;truth category;pairs", kTH1F, {axisTruthCategory});
    registry.add("LS/hPairTruthCategory", "like-sign pair truth category;truth category;pairs", kTH1F, {axisTruthCategory});
    registry.add("OS/hPairPtVsEta", "opposite-sign pair kinematics;p_{T}^{pair} (GeV/#it{c});#eta^{pair}", kTH2F, {axisPairPt, axisPairEta});
    registry.add("LS/hPairPtVsEta", "like-sign pair kinematics;p_{T}^{pair} (GeV/#it{c});#eta^{pair}", kTH2F, {axisPairPt, axisPairEta});

    addTruthMuonResonanceHistograms("Omega", axisMass, axisPairPt, axisPairEta);
    addTruthMuonResonanceHistograms("Phi", axisMass, axisPairPt, axisPairEta);
    addTruthMuonResonanceHistograms("Jpsi", axisMass, axisPairPt, axisPairEta);
  }

  void addPairHistograms(char const* hypothesis,
                         AxisSpec const& axisMass,
                         AxisSpec const& axisPairPt,
                         AxisSpec const& axisPairEta,
                         AxisSpec const& axisOpeningAngle,
                         AxisSpec const& axisTruthCategory)
  {
    const std::string base = hypothesis;
    registry.add("OS/hPairMass" + base, "opposite-sign pair mass;m (GeV/#it{c}^{2});pairs", kTH1F, {axisMass});
    registry.add("LS/hPairMass" + base, "like-sign pair mass;m (GeV/#it{c}^{2});pairs", kTH1F, {axisMass});
    registry.add("OS/hPairMassVsPt" + base, "opposite-sign pair mass vs pT;m (GeV/#it{c}^{2});p_{T}^{pair} (GeV/#it{c})", kTH2F, {axisMass, axisPairPt});
    registry.add("LS/hPairMassVsPt" + base, "like-sign pair mass vs pT;m (GeV/#it{c}^{2});p_{T}^{pair} (GeV/#it{c})", kTH2F, {axisMass, axisPairPt});
    registry.add("OS/hPairMassVsEta" + base, "opposite-sign pair mass vs eta;m (GeV/#it{c}^{2});#eta^{pair}", kTH2F, {axisMass, axisPairEta});
    registry.add("LS/hPairMassVsEta" + base, "like-sign pair mass vs eta;m (GeV/#it{c}^{2});#eta^{pair}", kTH2F, {axisMass, axisPairEta});
    registry.add("OS/hOpeningAngleVsMass" + base, "opposite-sign opening angle vs mass;m (GeV/#it{c}^{2});#alpha (rad)", kTH2F, {axisMass, axisOpeningAngle});
    registry.add("LS/hOpeningAngleVsMass" + base, "like-sign opening angle vs mass;m (GeV/#it{c}^{2});#alpha (rad)", kTH2F, {axisMass, axisOpeningAngle});
    registry.add("OS/hPairMassVsTruthCategory" + base, "opposite-sign mass vs truth category;m (GeV/#it{c}^{2});truth category", kTH2F, {axisMass, axisTruthCategory});
    registry.add("LS/hPairMassVsTruthCategory" + base, "like-sign mass vs truth category;m (GeV/#it{c}^{2});truth category", kTH2F, {axisMass, axisTruthCategory});
  }

  void addTruthMuonResonanceHistograms(char const* resonance,
                                       AxisSpec const& axisMass,
                                       AxisSpec const& axisPairPt,
                                       AxisSpec const& axisPairEta)
  {
    const std::string base = resonance;
    registry.add("TruthMuMu/" + base + "/hPairMassReco", "reconstructed mass for true resonance #rightarrow #mu#mu;m_{reco} (GeV/#it{c}^{2});pairs", kTH1F, {axisMass});
    registry.add("TruthMuMu/" + base + "/hPairMassTrue", "generated mass for true resonance #rightarrow #mu#mu;m_{true} (GeV/#it{c}^{2});pairs", kTH1F, {axisMass});
    registry.add("TruthMuMu/" + base + "/hPairMassRecoVsPt", "reconstructed mass vs pT for true resonance #rightarrow #mu#mu;m_{reco} (GeV/#it{c}^{2});p_{T}^{pair,reco} (GeV/#it{c})", kTH2F, {axisMass, axisPairPt});
    registry.add("TruthMuMu/" + base + "/hPairMassRecoVsEta", "reconstructed mass vs eta for true resonance #rightarrow #mu#mu;m_{reco} (GeV/#it{c}^{2});#eta^{pair,reco}", kTH2F, {axisMass, axisPairEta});
  }

  template <typename TAssoc>
  std::unordered_map<int64_t, int> countAssociationsPerTrack(TAssoc const& associations)
  {
    std::unordered_map<int64_t, int> nAssociationsPerTrack;
    nAssociationsPerTrack.reserve(associations.size());
    for (auto const& association : associations) {
      nAssociationsPerTrack[association.mfttrackId()]++;
    }
    return nAssociationsPerTrack;
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track) const
  {
    return track.pt() >= minTrackPt.value &&
           track.pt() <= maxTrackPt.value &&
           track.eta() >= minTrackEta.value &&
           track.eta() <= maxTrackEta.value &&
           track.nClusters() >= minNClusters.value &&
           track.chi2() <= maxTrackChi2.value &&
           std::isfinite(track.pt()) &&
           std::isfinite(track.eta()) &&
           std::isfinite(track.phi()) &&
           std::isfinite(track.tgl()) &&
           std::isfinite(track.signed1Pt());
  }

  static int getChargeSign(MFTTrackLabeled const& track)
  {
    if (track.signed1Pt() > 0.f) {
      return 1;
    }
    if (track.signed1Pt() < 0.f) {
      return -1;
    }
    return 0;
  }

  static FourVector makeFourVector(MFTTrackLabeled const& track, float mass)
  {
    const float pt = track.pt();
    const float px = pt * std::cos(track.phi());
    const float py = pt * std::sin(track.phi());
    const float pz = pt * track.tgl();
    const float p2 = px * px + py * py + pz * pz;
    return {px, py, pz, std::sqrt(p2 + mass * mass)};
  }

  static PairKinematics makePairKinematics(MFTTrackLabeled const& track1,
                                           MFTTrackLabeled const& track2,
                                           float mass1,
                                           float mass2)
  {
    const auto p1 = makeFourVector(track1, mass1);
    const auto p2 = makeFourVector(track2, mass2);

    const float px = p1.px + p2.px;
    const float py = p1.py + p2.py;
    const float pz = p1.pz + p2.pz;
    const float e = p1.e + p2.e;
    const float p2sum = px * px + py * py + pz * pz;
    const float mass2Pair = e * e - p2sum;
    const float pt = std::sqrt(px * px + py * py);
    const float p = std::sqrt(p2sum);
    const float eta = (p > std::abs(pz)) ? 0.5f * std::log((p + pz) / (p - pz)) : std::numeric_limits<float>::quiet_NaN();

    const float norm1 = std::sqrt(p1.px * p1.px + p1.py * p1.py + p1.pz * p1.pz);
    const float norm2 = std::sqrt(p2.px * p2.px + p2.py * p2.py + p2.pz * p2.pz);
    float openingAngle = std::numeric_limits<float>::quiet_NaN();
    if (norm1 > 0.f && norm2 > 0.f) {
      const float cosAngle = std::clamp((p1.px * p2.px + p1.py * p2.py + p1.pz * p2.pz) / (norm1 * norm2), -1.f, 1.f);
      openingAngle = std::acos(cosAngle);
    }

    return {std::sqrt(std::max(mass2Pair, 0.f)), pt, eta, openingAngle};
  }

  static PairKinematics makeTruePairKinematics(MFTTrackLabeled const& track1,
                                               MFTTrackLabeled const& track2,
                                               float mass1,
                                               float mass2)
  {
    if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
      return {};
    }

    auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
    auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
    const float px1 = mcParticle1.px();
    const float py1 = mcParticle1.py();
    const float pz1 = mcParticle1.pz();
    const float px2 = mcParticle2.px();
    const float py2 = mcParticle2.py();
    const float pz2 = mcParticle2.pz();
    const float e1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mass1 * mass1);
    const float e2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mass2 * mass2);

    const float px = px1 + px2;
    const float py = py1 + py2;
    const float pz = pz1 + pz2;
    const float e = e1 + e2;
    const float p2sum = px * px + py * py + pz * pz;
    const float mass2Pair = e * e - p2sum;
    const float pt = std::sqrt(px * px + py * py);
    const float p = std::sqrt(p2sum);
    const float eta = (p > std::abs(pz)) ? 0.5f * std::log((p + pz) / (p - pz)) : std::numeric_limits<float>::quiet_NaN();

    return {std::sqrt(std::max(mass2Pair, 0.f)), pt, eta, std::numeric_limits<float>::quiet_NaN()};
  }

  static int getAbsTruthPdg(MFTTrackLabeled const& track)
  {
    if (!track.has_mcParticle()) {
      return 0;
    }
    auto mcParticle = track.template mcParticle_as<aod::McParticles>();
    return std::abs(mcParticle.pdgCode());
  }

  static TruthCategory getTruthCategory(MFTTrackLabeled const& track1, MFTTrackLabeled const& track2)
  {
    const int pdg1 = getAbsTruthPdg(track1);
    const int pdg2 = getAbsTruthPdg(track2);
    if (pdg1 == 0 || pdg2 == 0) {
      return kUnknown;
    }
    if (pdg1 == 211 && pdg2 == 211) {
      return kPiPi;
    }
    if (pdg1 == 321 && pdg2 == 321) {
      return kKK;
    }
    if (pdg1 == 13 && pdg2 == 13) {
      return kMuMu;
    }
    if ((pdg1 == 211 && pdg2 == 321) || (pdg1 == 321 && pdg2 == 211)) {
      return kPiK;
    }
    return kOther;
  }

  static TruthMuonMother getTruthMuonMother(MFTTrackLabeled const& track)
  {
    if (!track.has_mcParticle()) {
      return {};
    }

    auto mcParticle = track.template mcParticle_as<aod::McParticles>();
    if (std::abs(mcParticle.pdgCode()) != 13 || !mcParticle.has_mothers()) {
      return {};
    }

    auto motherIds = mcParticle.mothersIds();
    if (motherIds.empty()) {
      return {};
    }

    auto mother = mcParticle.mothers_first_as<aod::McParticles>();
    return {true, std::abs(mother.pdgCode()), motherIds[0]};
  }

  static int getTrueMuonResonancePdg(MFTTrackLabeled const& track1, MFTTrackLabeled const& track2)
  {
    const auto mother1 = getTruthMuonMother(track1);
    const auto mother2 = getTruthMuonMother(track2);
    if (!mother1.valid || !mother2.valid || mother1.motherId != mother2.motherId || mother1.motherPdg != mother2.motherPdg) {
      return 0;
    }
    if (mother1.motherPdg == 223 || mother1.motherPdg == 333 || mother1.motherPdg == 443) {
      return mother1.motherPdg;
    }
    return 0;
  }

  void fillPairHistograms(bool isOppositeSign,
                          MassHypothesis hypothesis,
                          PairKinematics const& pair,
                          TruthCategory truthCategory)
  {
    if (!std::isfinite(pair.mass) || !std::isfinite(pair.pt) || !std::isfinite(pair.eta)) {
      return;
    }
    if (isOppositeSign) {
      if (hypothesis == kPiPiMass) {
        registry.fill(HIST("OS/hPairMassPiPi"), pair.mass);
        registry.fill(HIST("OS/hPairMassVsPtPiPi"), pair.mass, pair.pt);
        registry.fill(HIST("OS/hPairMassVsEtaPiPi"), pair.mass, pair.eta);
        registry.fill(HIST("OS/hPairMassVsTruthCategoryPiPi"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("OS/hOpeningAngleVsMassPiPi"), pair.mass, pair.openingAngle);
        }
      } else if (hypothesis == kKKMass) {
        registry.fill(HIST("OS/hPairMassKK"), pair.mass);
        registry.fill(HIST("OS/hPairMassVsPtKK"), pair.mass, pair.pt);
        registry.fill(HIST("OS/hPairMassVsEtaKK"), pair.mass, pair.eta);
        registry.fill(HIST("OS/hPairMassVsTruthCategoryKK"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("OS/hOpeningAngleVsMassKK"), pair.mass, pair.openingAngle);
        }
      } else {
        registry.fill(HIST("OS/hPairMassMuMu"), pair.mass);
        registry.fill(HIST("OS/hPairMassVsPtMuMu"), pair.mass, pair.pt);
        registry.fill(HIST("OS/hPairMassVsEtaMuMu"), pair.mass, pair.eta);
        registry.fill(HIST("OS/hPairMassVsTruthCategoryMuMu"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("OS/hOpeningAngleVsMassMuMu"), pair.mass, pair.openingAngle);
        }
      }
    } else {
      if (hypothesis == kPiPiMass) {
        registry.fill(HIST("LS/hPairMassPiPi"), pair.mass);
        registry.fill(HIST("LS/hPairMassVsPtPiPi"), pair.mass, pair.pt);
        registry.fill(HIST("LS/hPairMassVsEtaPiPi"), pair.mass, pair.eta);
        registry.fill(HIST("LS/hPairMassVsTruthCategoryPiPi"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("LS/hOpeningAngleVsMassPiPi"), pair.mass, pair.openingAngle);
        }
      } else if (hypothesis == kKKMass) {
        registry.fill(HIST("LS/hPairMassKK"), pair.mass);
        registry.fill(HIST("LS/hPairMassVsPtKK"), pair.mass, pair.pt);
        registry.fill(HIST("LS/hPairMassVsEtaKK"), pair.mass, pair.eta);
        registry.fill(HIST("LS/hPairMassVsTruthCategoryKK"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("LS/hOpeningAngleVsMassKK"), pair.mass, pair.openingAngle);
        }
      } else {
        registry.fill(HIST("LS/hPairMassMuMu"), pair.mass);
        registry.fill(HIST("LS/hPairMassVsPtMuMu"), pair.mass, pair.pt);
        registry.fill(HIST("LS/hPairMassVsEtaMuMu"), pair.mass, pair.eta);
        registry.fill(HIST("LS/hPairMassVsTruthCategoryMuMu"), pair.mass, static_cast<int>(truthCategory));
        if (std::isfinite(pair.openingAngle)) {
          registry.fill(HIST("LS/hOpeningAngleVsMassMuMu"), pair.mass, pair.openingAngle);
        }
      }
    }
  }

  void fillPairTruthAndKinematics(bool isOppositeSign, TruthCategory truthCategory, PairKinematics const& pair)
  {
    if (isOppositeSign) {
      registry.fill(HIST("OS/hPairTruthCategory"), static_cast<int>(truthCategory));
      if (std::isfinite(pair.pt) && std::isfinite(pair.eta)) {
        registry.fill(HIST("OS/hPairPtVsEta"), pair.pt, pair.eta);
      }
    } else {
      registry.fill(HIST("LS/hPairTruthCategory"), static_cast<int>(truthCategory));
      if (std::isfinite(pair.pt) && std::isfinite(pair.eta)) {
        registry.fill(HIST("LS/hPairPtVsEta"), pair.pt, pair.eta);
      }
    }
  }

  void fillTruthMuonResonanceHistograms(int motherPdg, PairKinematics const& recoPair, PairKinematics const& truePair)
  {
    if (!std::isfinite(recoPair.mass) || !std::isfinite(recoPair.pt) || !std::isfinite(recoPair.eta)) {
      return;
    }

    if (motherPdg == 223) {
      registry.fill(HIST("TruthMuMu/Omega/hPairMassReco"), recoPair.mass);
      registry.fill(HIST("TruthMuMu/Omega/hPairMassRecoVsPt"), recoPair.mass, recoPair.pt);
      registry.fill(HIST("TruthMuMu/Omega/hPairMassRecoVsEta"), recoPair.mass, recoPair.eta);
      if (std::isfinite(truePair.mass)) {
        registry.fill(HIST("TruthMuMu/Omega/hPairMassTrue"), truePair.mass);
      }
    } else if (motherPdg == 333) {
      registry.fill(HIST("TruthMuMu/Phi/hPairMassReco"), recoPair.mass);
      registry.fill(HIST("TruthMuMu/Phi/hPairMassRecoVsPt"), recoPair.mass, recoPair.pt);
      registry.fill(HIST("TruthMuMu/Phi/hPairMassRecoVsEta"), recoPair.mass, recoPair.eta);
      if (std::isfinite(truePair.mass)) {
        registry.fill(HIST("TruthMuMu/Phi/hPairMassTrue"), truePair.mass);
      }
    } else if (motherPdg == 443) {
      registry.fill(HIST("TruthMuMu/Jpsi/hPairMassReco"), recoPair.mass);
      registry.fill(HIST("TruthMuMu/Jpsi/hPairMassRecoVsPt"), recoPair.mass, recoPair.pt);
      registry.fill(HIST("TruthMuMu/Jpsi/hPairMassRecoVsEta"), recoPair.mass, recoPair.eta);
      if (std::isfinite(truePair.mass)) {
        registry.fill(HIST("TruthMuMu/Jpsi/hPairMassTrue"), truePair.mass);
      }
    }
  }

  void process(aod::Collisions const&,
               MFTTracksLabeled const&,
               aod::MFTTrackAssoc const& associations,
               aod::McParticles const&)
  {
    auto nAssociationsPerTrack = countAssociationsPerTrack(associations);

    std::unordered_map<int64_t, std::vector<MFTTrackLabeled>> tracksPerCollision;
    tracksPerCollision.reserve(associations.size());
    for (auto const& association : associations) {
      auto nAssoc = nAssociationsPerTrack[association.mfttrackId()];
      registry.fill(HIST("hAssociationMultiplicity"), nAssoc);
      if (nAssoc != 1) {
        continue;
      }

      auto track = association.mfttrack_as<MFTTracksLabeled>();
      if (!isSelectedTrack(track)) {
        continue;
      }

      registry.fill(HIST("hTrackPtEta"), track.pt(), track.eta());
      tracksPerCollision[association.collisionId()].push_back(track);
    }

    constexpr float massMuon = 0.1056583745f;
    constexpr float massPion = 0.13957039f;
    constexpr float massKaon = 0.493677f;

    for (auto& [collisionId, tracks] : tracksPerCollision) {
      (void)collisionId;
      registry.fill(HIST("hSelectedTracksPerCollision"), tracks.size());
      if (static_cast<int>(tracks.size()) > maxTracksPerCollision.value) {
        registry.fill(HIST("hSkippedCollisionsByMultiplicity"), tracks.size());
        continue;
      }

      for (std::size_t i = 0; i < tracks.size(); ++i) {
        for (std::size_t j = i + 1; j < tracks.size(); ++j) {
          auto const& track1 = tracks[i];
          auto const& track2 = tracks[j];
          const int sign1 = getChargeSign(track1);
          const int sign2 = getChargeSign(track2);
          if (sign1 == 0 || sign2 == 0) {
            continue;
          }

          const bool isOppositeSign = sign1 * sign2 < 0;
          const auto truthCategory = getTruthCategory(track1, track2);

          const auto pairPiPi = makePairKinematics(track1, track2, massPion, massPion);
          const auto pairKK = makePairKinematics(track1, track2, massKaon, massKaon);
          const auto pairMuMu = makePairKinematics(track1, track2, massMuon, massMuon);
          fillPairTruthAndKinematics(isOppositeSign, truthCategory, pairPiPi);
          fillPairHistograms(isOppositeSign, kPiPiMass, pairPiPi, truthCategory);
          fillPairHistograms(isOppositeSign, kKKMass, pairKK, truthCategory);
          fillPairHistograms(isOppositeSign, kMuMuMass, pairMuMu, truthCategory);

          const int trueMuonMotherPdg = getTrueMuonResonancePdg(track1, track2);
          if (trueMuonMotherPdg != 0) {
            const auto truePairMuMu = makeTruePairKinematics(track1, track2, massMuon, massMuon);
            fillTruthMuonResonanceHistograms(trueMuonMotherPdg, pairMuMu, truePairMuMu);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTTrackPairMass>(cfgc)};
}
