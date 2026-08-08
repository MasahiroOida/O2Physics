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

/// \file mftStandaloneTrackQa.cxx
/// \brief QA of MFT standalone tracks, with separate LTF and CA assessment.

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>

using namespace o2;
using namespace o2::framework;

struct MFTStandaloneTrackQa {
  using MFTTracksMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

  Service<o2::framework::O2DatabasePDG> pdg;

  struct : ConfigurableGroup {
    Configurable<float> etaMin{"cfgEtaMin", -3.6f, "Minimum track/truth eta (inclusive)"};
    Configurable<float> etaMax{"cfgEtaMax", -2.5f, "Maximum track/truth eta (inclusive)"};
    Configurable<float> ptMin{"cfgPtMin", 0.f, "Minimum track/truth pT in GeV/c (inclusive)"};
    Configurable<float> ptMax{"cfgPtMax", 1.e6f, "Maximum track/truth pT in GeV/c (inclusive)"};
    Configurable<int> minNClusters{"cfgMinNClusters", 5, "Minimum number of MFT clusters (inclusive)"};
    Configurable<float> maxChi2Ndf{"cfgMaxChi2Ndf", 1.e6f, "Maximum MFT track chi2/NDF (inclusive)"};
    Configurable<bool> physicalPrimaryOnly{"cfgPhysicalPrimaryOnly", true, "Use only physical-primary particles in the MC efficiency denominator"};
  } selection;

  ConfigurableAxis axisPt{"axisPt", {200, 0., 20.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisEta{"axisEta", {110, -4.5, -2.0}, "#it{#eta}"};
  ConfigurableAxis axisPhi{"axisPhi", {180, 0., o2::constants::math::TwoPI}, "#it{#varphi} (rad)"};
  ConfigurableAxis axisXY{"axisXY", {400, -20., 20.}, "position (cm)"};
  ConfigurableAxis axisZ{"axisZ", {200, -100., 0.}, "#it{z} (cm)"};
  ConfigurableAxis axisSigned1Pt{"axisSigned1Pt", {400, -20., 20.}, "#it{q}/#it{p}_{T} (#it{c}/GeV)"};
  ConfigurableAxis axisNClusters{"axisNClusters", {11, -0.5, 10.5}, "MFT clusters"};
  ConfigurableAxis axisClusterSize{"axisClusterSize", {64, -0.5, 63.5}, "cluster size"};
  ConfigurableAxis axisChi2Ndf{"axisChi2Ndf", {200, 0., 20.}, "#chi^{2}/NDF"};
  ConfigurableAxis axisTrackTime{"axisTrackTime", {400, -5000., 5000.}, "track time (ns)"};
  ConfigurableAxis axisTrackTimeRes{"axisTrackTimeRes", {200, 0., 1000.}, "track-time resolution (ns)"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {5000, -0.5, 4999.5}, "selected MFT tracks per collision"};
  ConfigurableAxis axisMultiplicityCoarse{"axisMultiplicityCoarse", {200, -0.5, 4999.5}, "selected MFT tracks per collision"};
  ConfigurableAxis axisPtRelativeResidual{"axisPtRelativeResidual", {400, -5., 5.}, "(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true}"};
  ConfigurableAxis axisEtaResidual{"axisEtaResidual", {400, -0.2, 0.2}, "#it{#eta}^{reco} - #it{#eta}^{true}"};
  ConfigurableAxis axisPhiResidual{"axisPhiResidual", {400, -0.2, 0.2}, "#it{#varphi}^{reco} - #it{#varphi}^{true} (rad)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct TrackCounts {
    int all = 0;
    int ltf = 0;
    int ca = 0;
  };

  enum RecoMask : uint8_t {
    kLTF = 0x1,
    kCA = 0x2
  };

  void addRecoHistograms()
  {
    registry.add("reco/all/pt", "MFT standalone tracks;#it{p}_{T} (GeV/#it{c});counts", HistType::kTH1F, {axisPt});
    registry.add("reco/all/eta", "MFT standalone tracks;#it{#eta};counts", HistType::kTH1F, {axisEta});
    registry.add("reco/all/phi", "MFT standalone tracks;#it{#varphi} (rad);counts", HistType::kTH1F, {axisPhi});
    registry.add("reco/all/z", "MFT standalone tracks;stored track #it{z} (cm);counts", HistType::kTH1F, {axisZ});
    registry.add("reco/all/signed1Pt", "MFT standalone tracks;#it{q}/#it{p}_{T} (#it{c}/GeV);counts", HistType::kTH1F, {axisSigned1Pt});
    registry.add("reco/all/sign", "MFT standalone tracks;charge sign;counts", HistType::kTH1F, {{3, -1.5, 1.5}});
    registry.add("reco/all/nClusters", "MFT standalone tracks;MFT clusters;counts", HistType::kTH1F, {axisNClusters});
    registry.add("reco/all/chi2Ndf", "MFT standalone tracks;#chi^{2}/NDF;counts", HistType::kTH1F, {axisChi2Ndf});
    registry.add("reco/all/xy", "MFT standalone tracks at stored track #it{z};#it{x} (cm);#it{y} (cm)", HistType::kTH2F, {axisXY, axisXY});
    registry.add("reco/all/etaPhi", "MFT standalone tracks;#it{#eta};#it{#varphi} (rad)", HistType::kTH2F, {axisEta, axisPhi});
    registry.add("reco/all/nClustersVsEta", "MFT standalone tracks;#it{#eta};MFT clusters", HistType::kTH2F, {axisEta, axisNClusters});
    registry.add("reco/all/chi2NdfVsPt", "MFT standalone tracks;#it{p}_{T} (GeV/#it{c});#chi^{2}/NDF", HistType::kTH2F, {axisPt, axisChi2Ndf});
    registry.add("reco/all/trackTime", "MFT standalone tracks;track time (ns);counts", HistType::kTH1F, {axisTrackTime});
    registry.add("reco/all/trackTimeRes", "MFT standalone tracks;track-time resolution (ns);counts", HistType::kTH1F, {axisTrackTimeRes});
    registry.add("reco/all/layerOccupancy", "MFT clusters used by standalone tracks;MFT layer;clusters", HistType::kTH1F, {{10, -0.5, 9.5}});
    registry.add("reco/all/clusterSizeVsLayer", "MFT cluster size by layer;MFT layer;cluster size", HistType::kTH2F, {{10, -0.5, 9.5}, axisClusterSize});
    registry.add("reco/all/layerPattern", "MFT standalone-track layer pattern;10-layer pattern;tracks", HistType::kTH1F, {{1024, -0.5, 1023.5}});
    registry.addClone("reco/all/", "reco/ltf/");
    registry.addClone("reco/all/", "reco/ca/");
  }

  void addCAContributionHistograms()
  {
    registry.add("caContribution/fractionVsPt", "CA contribution among selected tracks;#it{p}_{T}^{reco} (GeV/#it{c});CA fraction", HistType::kTProfile, {axisPt});
    registry.add("caContribution/fractionVsEta", "CA contribution among selected tracks;#it{#eta};CA fraction", HistType::kTProfile, {axisEta});
    registry.add("caContribution/fractionVsNClusters", "CA contribution among selected tracks;MFT clusters;CA fraction", HistType::kTProfile, {axisNClusters});
    registry.add("caContribution/fractionVsChi2Ndf", "CA contribution among selected tracks;#chi^{2}/NDF;CA fraction", HistType::kTProfile, {axisChi2Ndf});

    registry.add("collision/nTracksAll", "Selected MFT tracks per collision;MFT tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    registry.add("collision/nTracksLTF", "Selected LTF tracks per collision;LTF tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    registry.add("collision/nTracksCA", "Selected CA tracks per collision;CA tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    registry.add("collision/nCAVsNLTF", "CA and LTF multiplicities;LTF tracks;CA tracks", HistType::kTH2F, {axisMultiplicityCoarse, axisMultiplicityCoarse});
    registry.add("collision/caFraction", "CA contribution per collision;CA/(LTF+CA);collisions", HistType::kTH1F, {{110, -0.05, 1.05}});
    registry.add("collision/caFractionVsNTracks", "CA contribution vs selected multiplicity;selected MFT tracks;CA/(LTF+CA)", HistType::kTProfile, {axisMultiplicity});
  }

  void addMCHistograms()
  {
    registry.add("mc/labelStatusVsAlgorithm", "MC label status by tracking algorithm;tracking algorithm;MC label status", HistType::kTH2F, {{2, -0.5, 1.5}, {5, -0.5, 4.5}});

    registry.add("mc/ltf/nonCorrectLabelFractionVsPt", "Selected MFT tracks;#it{p}_{T}^{reco} (GeV/#it{c});non-correct label fraction", HistType::kTProfile, {axisPt});
    registry.add("mc/ltf/nonCorrectLabelFractionVsEta", "Selected MFT tracks;#it{#eta}^{reco};non-correct label fraction", HistType::kTProfile, {axisEta});
    registry.add("mc/ltf/nonCorrectLabelFractionVsNClusters", "Selected MFT tracks;MFT clusters;non-correct label fraction", HistType::kTProfile, {axisNClusters});
    registry.add("mc/ltf/ptRelativeResidualVsTruePt", "Selected truth-matched MFT tracks;#it{p}_{T}^{true} (GeV/#it{c});(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true}", HistType::kTH2F, {axisPt, axisPtRelativeResidual});
    registry.add("mc/ltf/etaResidualVsTrueEta", "Selected truth-matched MFT tracks;#it{#eta}^{true};#it{#eta}^{reco} - #it{#eta}^{true}", HistType::kTH2F, {axisEta, axisEtaResidual});
    registry.add("mc/ltf/phiResidualVsTrueEta", "Selected truth-matched MFT tracks;#it{#eta}^{true};#it{#varphi}^{reco} - #it{#varphi}^{true} (rad)", HistType::kTH2F, {axisEta, axisPhiResidual});
    registry.add("mc/ltf/chargeCorrectVsTruePt", "Selected truth-matched MFT tracks;#it{p}_{T}^{true} (GeV/#it{c});correct charge fraction", HistType::kTProfile, {axisPt});
    registry.addClone("mc/ltf/", "mc/ca/");

    registry.add("mc/truth/denominator/pt", "Selected charged truth particles;#it{p}_{T}^{true} (GeV/#it{c});particles", HistType::kTH1F, {axisPt});
    registry.add("mc/truth/denominator/eta", "Selected charged truth particles;#it{#eta}^{true};particles", HistType::kTH1F, {axisEta});
    registry.add("mc/truth/denominator/ptEta", "Selected charged truth particles;#it{p}_{T}^{true} (GeV/#it{c});#it{#eta}^{true}", HistType::kTH2F, {axisPt, axisEta});
    registry.addClone("mc/truth/denominator/", "mc/truth/inclusive/");
    registry.addClone("mc/truth/denominator/", "mc/truth/ltf/");
    registry.addClone("mc/truth/denominator/", "mc/truth/ca/");
    registry.addClone("mc/truth/denominator/", "mc/truth/ltfOnly/");
    registry.addClone("mc/truth/denominator/", "mc/truth/caOnly/");
    registry.addClone("mc/truth/denominator/", "mc/truth/both/");
    registry.add("mc/truth/recoveryCategory", "Selected truth-particle recovery category;recovery category;particles", HistType::kTH1F, {{4, -0.5, 3.5}});
    registry.add("mc/truth/nRecoTracks", "Selected reconstructed tracks per truth particle;reconstructed tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
    registry.add("mc/truth/nRecoTracksLTF", "Selected LTF tracks per truth particle;LTF tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
    registry.add("mc/truth/nRecoTracksCA", "Selected CA tracks per truth particle;CA tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
  }

  void init(InitContext const&)
  {
    if (doprocessData == doprocessMC) {
      LOGF(fatal, "Enable exactly one of processData and processMC");
    }
    if (selection.etaMin.value > selection.etaMax.value || selection.ptMin.value > selection.ptMax.value) {
      LOGF(fatal, "Invalid MFT standalone QA selection range");
    }
    if (selection.minNClusters.value < 0 || selection.maxChi2Ndf.value < 0.f) {
      LOGF(fatal, "MFT cluster and chi2/NDF selections must be non-negative");
    }

    addRecoHistograms();
    addCAContributionHistograms();
    if (doprocessMC) {
      addMCHistograms();
    }

    auto selectionCounter = registry.add<TH1>("reco/selectionCounter", "MFT standalone track selection;;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    selectionCounter->GetXaxis()->SetBinLabel(1, "all");
    selectionCounter->GetXaxis()->SetBinLabel(2, "selected");

    auto algorithm = registry.add<TH1>("reco/selectedAlgorithm", "Selected MFT standalone tracks;tracking algorithm;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    algorithm->GetXaxis()->SetBinLabel(1, "LTF");
    algorithm->GetXaxis()->SetBinLabel(2, "CA");

    auto association = registry.add<TH1>("reco/associationStatus", "MFT standalone track collision association;association;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    association->GetXaxis()->SetBinLabel(1, "unassociated");
    association->GetXaxis()->SetBinLabel(2, "associated");

    if (doprocessMC) {
      auto labelStatus = registry.get<TH2>(HIST("mc/labelStatusVsAlgorithm"));
      labelStatus->GetXaxis()->SetBinLabel(1, "LTF");
      labelStatus->GetXaxis()->SetBinLabel(2, "CA");
      labelStatus->GetYaxis()->SetBinLabel(1, "correct");
      labelStatus->GetYaxis()->SetBinLabel(2, "fake");
      labelStatus->GetYaxis()->SetBinLabel(3, "noise");
      labelStatus->GetYaxis()->SetBinLabel(4, "fake+noise");
      labelStatus->GetYaxis()->SetBinLabel(5, "no label");

      auto recovery = registry.get<TH1>(HIST("mc/truth/recoveryCategory"));
      recovery->GetXaxis()->SetBinLabel(1, "missed");
      recovery->GetXaxis()->SetBinLabel(2, "LTF only");
      recovery->GetXaxis()->SetBinLabel(3, "CA only");
      recovery->GetXaxis()->SetBinLabel(4, "LTF and CA");
    }
  }

  template <typename TTrack>
  float chi2Ndf(TTrack const& track) const
  {
    const int ndf = 2 * static_cast<int>(track.nClusters()) - 5;
    if (ndf <= 0 || !std::isfinite(track.chi2())) {
      return -1.f;
    }
    return track.chi2() / static_cast<float>(ndf);
  }

  float wrapPhi(float phi) const
  {
    float wrapped = std::fmod(phi, static_cast<float>(o2::constants::math::TwoPI));
    if (wrapped < 0.f) {
      wrapped += static_cast<float>(o2::constants::math::TwoPI);
    }
    return wrapped;
  }

  float deltaPhi(float phiReco, float phiTrue) const
  {
    return std::atan2(std::sin(phiReco - phiTrue), std::cos(phiReco - phiTrue));
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track) const
  {
    const float normalizedChi2 = chi2Ndf(track);
    return std::isfinite(track.pt()) && std::isfinite(track.eta()) && normalizedChi2 >= 0.f &&
           track.pt() >= selection.ptMin.value && track.pt() <= selection.ptMax.value &&
           track.eta() >= selection.etaMin.value && track.eta() <= selection.etaMax.value &&
           track.nClusters() >= selection.minNClusters.value && normalizedChi2 <= selection.maxChi2Ndf.value;
  }

  template <typename TParticle>
  bool isSelectedTruthParticle(TParticle const& particle) const
  {
    if (!std::isfinite(particle.pt()) || !std::isfinite(particle.eta()) ||
        particle.pt() < selection.ptMin.value || particle.pt() > selection.ptMax.value ||
        particle.eta() < selection.etaMin.value || particle.eta() > selection.etaMax.value) {
      return false;
    }
    if (selection.physicalPrimaryOnly.value && !particle.isPhysicalPrimary()) {
      return false;
    }
    const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
    return pdgParticle != nullptr && std::abs(pdgParticle->Charge()) >= 3.;
  }

  template <typename TTrack>
  void fillAllTrackHistograms(TTrack const& track, float normalizedChi2, float phi)
  {
    registry.fill(HIST("reco/all/pt"), track.pt());
    registry.fill(HIST("reco/all/eta"), track.eta());
    registry.fill(HIST("reco/all/phi"), phi);
    registry.fill(HIST("reco/all/z"), track.z());
    registry.fill(HIST("reco/all/signed1Pt"), track.signed1Pt());
    registry.fill(HIST("reco/all/sign"), track.sign());
    registry.fill(HIST("reco/all/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      registry.fill(HIST("reco/all/chi2Ndf"), normalizedChi2);
      registry.fill(HIST("reco/all/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    registry.fill(HIST("reco/all/xy"), track.x(), track.y());
    registry.fill(HIST("reco/all/etaPhi"), track.eta(), phi);
    registry.fill(HIST("reco/all/nClustersVsEta"), track.eta(), track.nClusters());
    registry.fill(HIST("reco/all/trackTime"), track.trackTime());
    registry.fill(HIST("reco/all/trackTimeRes"), track.trackTimeRes());
  }

  template <typename TTrack>
  void fillAlgorithmTrackHistograms(TTrack const& track, float normalizedChi2, float phi)
  {
    if (track.isCA()) {
      registry.fill(HIST("reco/ca/pt"), track.pt());
      registry.fill(HIST("reco/ca/eta"), track.eta());
      registry.fill(HIST("reco/ca/phi"), phi);
      registry.fill(HIST("reco/ca/z"), track.z());
      registry.fill(HIST("reco/ca/signed1Pt"), track.signed1Pt());
      registry.fill(HIST("reco/ca/sign"), track.sign());
      registry.fill(HIST("reco/ca/nClusters"), track.nClusters());
      if (normalizedChi2 >= 0.f) {
        registry.fill(HIST("reco/ca/chi2Ndf"), normalizedChi2);
        registry.fill(HIST("reco/ca/chi2NdfVsPt"), track.pt(), normalizedChi2);
      }
      registry.fill(HIST("reco/ca/xy"), track.x(), track.y());
      registry.fill(HIST("reco/ca/etaPhi"), track.eta(), phi);
      registry.fill(HIST("reco/ca/nClustersVsEta"), track.eta(), track.nClusters());
      registry.fill(HIST("reco/ca/trackTime"), track.trackTime());
      registry.fill(HIST("reco/ca/trackTimeRes"), track.trackTimeRes());
      return;
    }

    registry.fill(HIST("reco/ltf/pt"), track.pt());
    registry.fill(HIST("reco/ltf/eta"), track.eta());
    registry.fill(HIST("reco/ltf/phi"), phi);
    registry.fill(HIST("reco/ltf/z"), track.z());
    registry.fill(HIST("reco/ltf/signed1Pt"), track.signed1Pt());
    registry.fill(HIST("reco/ltf/sign"), track.sign());
    registry.fill(HIST("reco/ltf/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      registry.fill(HIST("reco/ltf/chi2Ndf"), normalizedChi2);
      registry.fill(HIST("reco/ltf/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    registry.fill(HIST("reco/ltf/xy"), track.x(), track.y());
    registry.fill(HIST("reco/ltf/etaPhi"), track.eta(), phi);
    registry.fill(HIST("reco/ltf/nClustersVsEta"), track.eta(), track.nClusters());
    registry.fill(HIST("reco/ltf/trackTime"), track.trackTime());
    registry.fill(HIST("reco/ltf/trackTimeRes"), track.trackTimeRes());
  }

  template <typename TTrack>
  void fillClusterLayerHistograms(TTrack const& track)
  {
    const uint64_t clusterSizesAndFlags = track.mftClusterSizesAndTrackFlags();
    uint16_t layerPattern = 0;
    for (int layer = 0; layer < 10; ++layer) {
      const int clusterSize = static_cast<int>((clusterSizesAndFlags >> (6 * layer)) & 0x3fULL);
      if (clusterSize == 0) {
        continue;
      }
      layerPattern |= static_cast<uint16_t>(1U << layer);
      registry.fill(HIST("reco/all/layerOccupancy"), layer);
      registry.fill(HIST("reco/all/clusterSizeVsLayer"), layer, clusterSize);
      if (track.isCA()) {
        registry.fill(HIST("reco/ca/layerOccupancy"), layer);
        registry.fill(HIST("reco/ca/clusterSizeVsLayer"), layer, clusterSize);
      } else {
        registry.fill(HIST("reco/ltf/layerOccupancy"), layer);
        registry.fill(HIST("reco/ltf/clusterSizeVsLayer"), layer, clusterSize);
      }
    }
    registry.fill(HIST("reco/all/layerPattern"), layerPattern);
    if (track.isCA()) {
      registry.fill(HIST("reco/ca/layerPattern"), layerPattern);
    } else {
      registry.fill(HIST("reco/ltf/layerPattern"), layerPattern);
    }
  }

  template <typename TTracks>
  void fillRecoQA(aod::Collisions const& collisions, TTracks const& tracks)
  {
    std::unordered_map<int64_t, TrackCounts> countsPerCollision;

    for (const auto& track : tracks) {
      registry.fill(HIST("reco/selectionCounter"), 0.);
      registry.fill(HIST("reco/associationStatus"), track.collisionId() < 0 ? 0. : 1.);

      const float normalizedChi2 = chi2Ndf(track);
      const float phi = wrapPhi(track.phi());
      fillAllTrackHistograms(track, normalizedChi2, phi);
      fillAlgorithmTrackHistograms(track, normalizedChi2, phi);
      fillClusterLayerHistograms(track);

      if (!isSelectedTrack(track)) {
        continue;
      }
      registry.fill(HIST("reco/selectionCounter"), 1.);
      registry.fill(HIST("reco/selectedAlgorithm"), track.isCA() ? 1. : 0.);
      registry.fill(HIST("caContribution/fractionVsPt"), track.pt(), track.isCA());
      registry.fill(HIST("caContribution/fractionVsEta"), track.eta(), track.isCA());
      registry.fill(HIST("caContribution/fractionVsNClusters"), track.nClusters(), track.isCA());
      registry.fill(HIST("caContribution/fractionVsChi2Ndf"), normalizedChi2, track.isCA());

      if (track.collisionId() >= 0) {
        auto& counts = countsPerCollision[track.collisionId()];
        ++counts.all;
        if (track.isCA()) {
          ++counts.ca;
        } else {
          ++counts.ltf;
        }
      }
    }

    for (const auto& collision : collisions) {
      const auto found = countsPerCollision.find(collision.globalIndex());
      const TrackCounts counts = found == countsPerCollision.end() ? TrackCounts{} : found->second;
      registry.fill(HIST("collision/nTracksAll"), counts.all);
      registry.fill(HIST("collision/nTracksLTF"), counts.ltf);
      registry.fill(HIST("collision/nTracksCA"), counts.ca);
      registry.fill(HIST("collision/nCAVsNLTF"), counts.ltf, counts.ca);
      if (counts.all > 0) {
        const float caFraction = static_cast<float>(counts.ca) / static_cast<float>(counts.all);
        registry.fill(HIST("collision/caFraction"), caFraction);
        registry.fill(HIST("collision/caFractionVsNTracks"), counts.all, caFraction);
      }
    }
  }

  template <typename TParticle>
  void fillTruthKinematics(char category, TParticle const& particle)
  {
    switch (category) {
      case 'd':
        registry.fill(HIST("mc/truth/denominator/pt"), particle.pt());
        registry.fill(HIST("mc/truth/denominator/eta"), particle.eta());
        registry.fill(HIST("mc/truth/denominator/ptEta"), particle.pt(), particle.eta());
        break;
      case 'i':
        registry.fill(HIST("mc/truth/inclusive/pt"), particle.pt());
        registry.fill(HIST("mc/truth/inclusive/eta"), particle.eta());
        registry.fill(HIST("mc/truth/inclusive/ptEta"), particle.pt(), particle.eta());
        break;
      case 'l':
        registry.fill(HIST("mc/truth/ltf/pt"), particle.pt());
        registry.fill(HIST("mc/truth/ltf/eta"), particle.eta());
        registry.fill(HIST("mc/truth/ltf/ptEta"), particle.pt(), particle.eta());
        break;
      case 'c':
        registry.fill(HIST("mc/truth/ca/pt"), particle.pt());
        registry.fill(HIST("mc/truth/ca/eta"), particle.eta());
        registry.fill(HIST("mc/truth/ca/ptEta"), particle.pt(), particle.eta());
        break;
      case 'L':
        registry.fill(HIST("mc/truth/ltfOnly/pt"), particle.pt());
        registry.fill(HIST("mc/truth/ltfOnly/eta"), particle.eta());
        registry.fill(HIST("mc/truth/ltfOnly/ptEta"), particle.pt(), particle.eta());
        break;
      case 'C':
        registry.fill(HIST("mc/truth/caOnly/pt"), particle.pt());
        registry.fill(HIST("mc/truth/caOnly/eta"), particle.eta());
        registry.fill(HIST("mc/truth/caOnly/ptEta"), particle.pt(), particle.eta());
        break;
      case 'b':
        registry.fill(HIST("mc/truth/both/pt"), particle.pt());
        registry.fill(HIST("mc/truth/both/eta"), particle.eta());
        registry.fill(HIST("mc/truth/both/ptEta"), particle.pt(), particle.eta());
        break;
    }
  }

  template <typename TTrack, typename TParticle>
  void fillResolution(TTrack const& track, TParticle const& particle)
  {
    const float phiResidual = deltaPhi(track.phi(), particle.phi());
    const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
    const bool chargeCorrect = pdgParticle != nullptr &&
                               ((pdgParticle->Charge() > 0. && track.sign() > 0) || (pdgParticle->Charge() < 0. && track.sign() < 0));

    if (track.isCA()) {
      if (particle.pt() > 0.f) {
        registry.fill(HIST("mc/ca/ptRelativeResidualVsTruePt"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
      }
      registry.fill(HIST("mc/ca/etaResidualVsTrueEta"), particle.eta(), track.eta() - particle.eta());
      registry.fill(HIST("mc/ca/phiResidualVsTrueEta"), particle.eta(), phiResidual);
      registry.fill(HIST("mc/ca/chargeCorrectVsTruePt"), particle.pt(), chargeCorrect);
      return;
    }

    if (particle.pt() > 0.f) {
      registry.fill(HIST("mc/ltf/ptRelativeResidualVsTruePt"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
    }
    registry.fill(HIST("mc/ltf/etaResidualVsTrueEta"), particle.eta(), track.eta() - particle.eta());
    registry.fill(HIST("mc/ltf/phiResidualVsTrueEta"), particle.eta(), phiResidual);
    registry.fill(HIST("mc/ltf/chargeCorrectVsTruePt"), particle.pt(), chargeCorrect);
  }

  template <typename TTrack>
  int labelStatus(TTrack const& track) const
  {
    constexpr uint8_t noiseBit = 0x1 << 6;
    constexpr uint8_t fakeBit = 0x1 << 7;
    const bool isNoise = track.mcMask() & noiseBit;
    const bool isFake = track.mcMask() & fakeBit;
    if (isNoise && isFake) {
      return 3;
    }
    if (isFake) {
      return 1;
    }
    if (isNoise) {
      return 2;
    }
    return track.has_mcParticle() ? 0 : 4;
  }

  void fillMCQA(MFTTracksMC const& tracks, aod::McParticles const& particles)
  {
    std::unordered_map<int64_t, uint8_t> recoMaskByParticle;
    std::unordered_map<int64_t, std::array<int, 2>> recoCountsByParticle;

    for (const auto& track : tracks) {
      if (!isSelectedTrack(track)) {
        continue;
      }

      const int status = labelStatus(track);
      const bool hasNonCorrectLabel = status != 0;
      registry.fill(HIST("mc/labelStatusVsAlgorithm"), track.isCA() ? 1. : 0., status);
      if (track.isCA()) {
        registry.fill(HIST("mc/ca/nonCorrectLabelFractionVsPt"), track.pt(), hasNonCorrectLabel);
        registry.fill(HIST("mc/ca/nonCorrectLabelFractionVsEta"), track.eta(), hasNonCorrectLabel);
        registry.fill(HIST("mc/ca/nonCorrectLabelFractionVsNClusters"), track.nClusters(), hasNonCorrectLabel);
      } else {
        registry.fill(HIST("mc/ltf/nonCorrectLabelFractionVsPt"), track.pt(), hasNonCorrectLabel);
        registry.fill(HIST("mc/ltf/nonCorrectLabelFractionVsEta"), track.eta(), hasNonCorrectLabel);
        registry.fill(HIST("mc/ltf/nonCorrectLabelFractionVsNClusters"), track.nClusters(), hasNonCorrectLabel);
      }

      if (status != 0) {
        continue;
      }
      const auto particle = track.template mcParticle_as<aod::McParticles>();
      if (!isSelectedTruthParticle(particle)) {
        continue;
      }

      auto& recoMask = recoMaskByParticle[particle.globalIndex()];
      auto& recoCounts = recoCountsByParticle[particle.globalIndex()];
      if (track.isCA()) {
        recoMask |= kCA;
        ++recoCounts[1];
      } else {
        recoMask |= kLTF;
        ++recoCounts[0];
      }
      fillResolution(track, particle);
    }

    for (const auto& particle : particles) {
      if (!isSelectedTruthParticle(particle)) {
        continue;
      }

      fillTruthKinematics('d', particle);
      const auto maskFound = recoMaskByParticle.find(particle.globalIndex());
      const uint8_t recoMask = maskFound == recoMaskByParticle.end() ? 0 : maskFound->second;
      const auto countsFound = recoCountsByParticle.find(particle.globalIndex());
      const std::array<int, 2> recoCounts = countsFound == recoCountsByParticle.end() ? std::array<int, 2>{0, 0} : countsFound->second;
      registry.fill(HIST("mc/truth/nRecoTracks"), recoCounts[0] + recoCounts[1]);
      registry.fill(HIST("mc/truth/nRecoTracksLTF"), recoCounts[0]);
      registry.fill(HIST("mc/truth/nRecoTracksCA"), recoCounts[1]);

      if (recoMask == 0) {
        registry.fill(HIST("mc/truth/recoveryCategory"), 0.);
        continue;
      }

      fillTruthKinematics('i', particle);
      if (recoMask & kLTF) {
        fillTruthKinematics('l', particle);
      }
      if (recoMask & kCA) {
        fillTruthKinematics('c', particle);
      }
      if (recoMask == kLTF) {
        fillTruthKinematics('L', particle);
        registry.fill(HIST("mc/truth/recoveryCategory"), 1.);
      } else if (recoMask == kCA) {
        fillTruthKinematics('C', particle);
        registry.fill(HIST("mc/truth/recoveryCategory"), 2.);
      } else {
        fillTruthKinematics('b', particle);
        registry.fill(HIST("mc/truth/recoveryCategory"), 3.);
      }
    }
  }

  void processData(aod::Collisions const& collisions, aod::MFTTracks const& tracks)
  {
    fillRecoQA(collisions, tracks);
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processData, "Process data MFT standalone-track QA", true);

  void processMC(aod::Collisions const& collisions, MFTTracksMC const& tracks, aod::McParticles const& particles)
  {
    fillRecoQA(collisions, tracks);
    fillMCQA(tracks, particles);
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processMC, "Process MC MFT standalone-track QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTStandaloneTrackQa>(cfgc)};
}
