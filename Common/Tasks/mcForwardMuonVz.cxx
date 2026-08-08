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

/// \file mcForwardMuonVz.cxx
/// \brief QA task for the production-z distribution of generated forward muons.

#include <Framework/AnalysisDataModel.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TMCProcess.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;

using LabeledFwdTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using LabeledMftTracks = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

struct McForwardMuonVz {
  Configurable<float> cfgEtaMin{"cfgEtaMin", -3.6f, "Minimum generated-muon pseudorapidity"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", -2.5f, "Maximum generated-muon pseudorapidity"};
  ConfigurableAxis cfgVzAxis{"cfgVzAxis", {400, -800.f, 800.f}, "Generated-muon z vertex axis in cm"};
  ConfigurableAxis cfgTruePtAxis{"cfgTruePtAxis", {100, 0.f, 10.f}, "True-pT axis in GeV/c"};
  ConfigurableAxis cfgPiKDeltaPtAxis{"cfgPiKDeltaPtAxis", {400, -20.f, 20.f}, "#Delta pT axis for #pi/K #rightarrow #mu"};
  ConfigurableAxis cfgPiKDeltaEtaAxis{"cfgPiKDeltaEtaAxis", {400, -4.f, 4.f}, "#Delta #eta axis for #pi/K #rightarrow #mu"};
  ConfigurableAxis cfgPiKDeltaPhiAxis{"cfgPiKDeltaPhiAxis", {360, -3.1415927f, 3.1415927f}, "#Delta #phi axis for #pi/K #rightarrow #mu"};
  ConfigurableAxis cfgPiKDeltaVzAxis{"cfgPiKDeltaVzAxis", {400, -800.f, 800.f}, "#Delta z axis for #pi/K #rightarrow #mu in cm"};
  ConfigurableAxis cfgPiKDecayLengthAxis{"cfgPiKDecayLengthAxis", {500, 0.f, 50000.f}, "3D decay-length axis for #pi/K #rightarrow #mu in cm"};
  ConfigurableAxis cfgProcessAxis{"cfgProcessAxis", {51, -0.5f, 50.5f}, "TMCProcess code axis"};
  ConfigurableAxis cfgPdgAxis{"cfgPdgAxis", {10001, -5000.5f, 5000.5f}, "PDG code axis"};

  HistogramRegistry registry{
    "registry",
    {
      // Generated forward-muon production-z categories.
      {"forwardMuonVz", "Generated #mu^{#pm} production z for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzGenerator", "Generator #mu^{#pm} production z for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzNaturalDecay", "#mu^{#pm} from natural particle decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzPiKDecay", "#mu^{#pm} from natural #pi/K decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzMaterial", "#mu^{#pm} from material interaction for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzOther", "Other #mu^{#pm} production z for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"muonVzVsMuonProcess", "Generated #mu^{#pm} production z by muon process;z_{#mu} (cm);muon TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"muonVzVsMotherProcess", "Generated #mu^{#pm} production z by direct-mother process;z_{#mu} (cm);mother TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"motherVzVsMotherProcess", "Direct-mother production z by mother process for generated #mu^{#pm};z_{mother} (cm);mother TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},

      // Direct pi/K -> mu decays, split by the pi/K mother origin.
      {"forwardMuonVzPiKFromGenerator", "#mu^{#pm} from generator #pi/K decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzPiKFromMaterial", "#mu^{#pm} from material-origin #pi/K decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzPiKFromDecay", "#mu^{#pm} from decay-origin #pi/K decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"forwardMuonVzPiKFromOther", "#mu^{#pm} from other-origin #pi/K decay for -3.6 < #eta < -2.5;z_{gen} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherVzFromGenerator", "Generator #pi/K mother production z for #pi/K #rightarrow #mu;z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherVzFromMaterial", "Material-origin #pi/K mother production z for #pi/K #rightarrow #mu;z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherVzFromDecay", "Decay-origin #pi/K mother production z for #pi/K #rightarrow #mu;z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherVzFromOther", "Other-origin #pi/K mother production z for #pi/K #rightarrow #mu;z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},

      // Direct pi/K -> mu decay kinematic differences: muon - mother.
      {"piKToMuonDeltaPt", "#mu - #pi/K mother for #pi/K #rightarrow #mu;#Delta p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {cfgPiKDeltaPtAxis}}},
      {"piKToMuonDeltaEta", "#mu - #pi/K mother for #pi/K #rightarrow #mu;#Delta #eta;counts", {HistType::kTH1F, {cfgPiKDeltaEtaAxis}}},
      {"piKToMuonDeltaPhi", "#mu - #pi/K mother for #pi/K #rightarrow #mu;#Delta #phi;counts", {HistType::kTH1F, {cfgPiKDeltaPhiAxis}}},
      {"piKToMuonDeltaVz", "#mu - #pi/K mother for #pi/K #rightarrow #mu;#Delta z = z_{#mu} - z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgPiKDeltaVzAxis}}},
      {"piKToMuonDecayLength", "#pi/K #rightarrow #mu 3D decay length;L_{decay} (cm);counts", {HistType::kTH1F, {cfgPiKDecayLengthAxis}}},
      {"piKToMuonDecayLengthVsMuonVz", "#pi/K #rightarrow #mu 3D decay length vs #mu production z;z_{#mu} (cm);L_{decay} (cm);counts", {HistType::kTH2F, {cfgVzAxis, cfgPiKDecayLengthAxis}}},
      {"piToMuonDeltaPt", "#mu - #pi mother for #pi #rightarrow #mu;#Delta p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {cfgPiKDeltaPtAxis}}},
      {"piToMuonDeltaEta", "#mu - #pi mother for #pi #rightarrow #mu;#Delta #eta;counts", {HistType::kTH1F, {cfgPiKDeltaEtaAxis}}},
      {"piToMuonDeltaPhi", "#mu - #pi mother for #pi #rightarrow #mu;#Delta #phi;counts", {HistType::kTH1F, {cfgPiKDeltaPhiAxis}}},
      {"piToMuonDeltaVz", "#mu - #pi mother for #pi #rightarrow #mu;#Delta z = z_{#mu} - z_{#pi} (cm);counts", {HistType::kTH1F, {cfgPiKDeltaVzAxis}}},
      {"piToMuonDecayLength", "#pi #rightarrow #mu 3D decay length;L_{decay} (cm);counts", {HistType::kTH1F, {cfgPiKDecayLengthAxis}}},
      {"piToMuonDecayLengthVsMuonVz", "#pi #rightarrow #mu 3D decay length vs #mu production z;z_{#mu} (cm);L_{decay} (cm);counts", {HistType::kTH2F, {cfgVzAxis, cfgPiKDecayLengthAxis}}},
      {"kToMuonDeltaPt", "#mu - K mother for K #rightarrow #mu;#Delta p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {cfgPiKDeltaPtAxis}}},
      {"kToMuonDeltaEta", "#mu - K mother for K #rightarrow #mu;#Delta #eta;counts", {HistType::kTH1F, {cfgPiKDeltaEtaAxis}}},
      {"kToMuonDeltaPhi", "#mu - K mother for K #rightarrow #mu;#Delta #phi;counts", {HistType::kTH1F, {cfgPiKDeltaPhiAxis}}},
      {"kToMuonDeltaVz", "#mu - K mother for K #rightarrow #mu;#Delta z = z_{#mu} - z_{K} (cm);counts", {HistType::kTH1F, {cfgPiKDeltaVzAxis}}},
      {"kToMuonDecayLength", "K #rightarrow #mu 3D decay length;L_{decay} (cm);counts", {HistType::kTH1F, {cfgPiKDecayLengthAxis}}},
      {"kToMuonDecayLengthVsMuonVz", "K #rightarrow #mu 3D decay length vs #mu production z;z_{#mu} (cm);L_{decay} (cm);counts", {HistType::kTH2F, {cfgVzAxis, cfgPiKDecayLengthAxis}}},

      // All direct daughters of pi/K mothers.
      {"piKMotherDaughterVz", "Production z of particles with direct #pi/K mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherDaughterVzVsProcess", "Production z of particles with direct #pi/K mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"piKMotherNonMuonDaughterVz", "Production z of non-#mu particles with direct #pi/K mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piKMotherNonMuonDaughterVzVsProcess", "Production z of non-#mu particles with direct #pi/K mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"piKMotherNonMuonDaughterVzVsPdg", "Production z of non-#mu particles with direct #pi/K mother by child PDG;z_{child} (cm);child PDG;counts", {HistType::kTH2F, {cfgVzAxis, cfgPdgAxis}}},
      {"piMotherDaughterVz", "Production z of particles with direct #pi mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piMotherDaughterVzVsProcess", "Production z of particles with direct #pi mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"piMotherNonMuonDaughterVz", "Production z of non-#mu particles with direct #pi mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"piMotherNonMuonDaughterVzVsProcess", "Production z of non-#mu particles with direct #pi mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"piMotherNonMuonDaughterVzVsPdg", "Production z of non-#mu particles with direct #pi mother by child PDG;z_{child} (cm);child PDG;counts", {HistType::kTH2F, {cfgVzAxis, cfgPdgAxis}}},
      {"kMotherDaughterVz", "Production z of particles with direct K mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"kMotherDaughterVzVsProcess", "Production z of particles with direct K mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"kMotherNonMuonDaughterVz", "Production z of non-#mu particles with direct K mother;z_{child} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"kMotherNonMuonDaughterVzVsProcess", "Production z of non-#mu particles with direct K mother by child process;z_{child} (cm);child TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"kMotherNonMuonDaughterVzVsPdg", "Production z of non-#mu particles with direct K mother by child PDG;z_{child} (cm);child PDG;counts", {HistType::kTH2F, {cfgVzAxis, cfgPdgAxis}}},

      // Reconstructed-track truth counters and FwdTrack purity.
      {"muonTrackCounts", "MC-labeled muon counts for -3.6 < #eta < -2.5;;counts", {HistType::kTH1F, {{5, 0.5, 5.5}}}},
      {"fwdTrackTruePtAll", "All MC-labeled FwdTracks for -3.6 < #eta_{true} < -2.5;p_{T,true} (GeV/#it{c});counts", {HistType::kTH1F, {cfgTruePtAxis}}},
      {"fwdTrackTruePtMuon", "Muon MC-labeled FwdTracks for -3.6 < #eta_{true} < -2.5;p_{T,true} (GeV/#it{c});counts", {HistType::kTH1F, {cfgTruePtAxis}}},
      {"fwdTrackMuonFractionVsTruePt", "Muon fraction of MC-labeled FwdTracks for -3.6 < #eta_{true} < -2.5;p_{T,true} (GeV/#it{c});muon fraction", {HistType::kTProfile, {cfgTruePtAxis}}},

      // MFT-track MC-label production-z diagnostics.
      {"mftMuonProductionVz", "MFT-labeled #mu^{#pm} production z;z_{#mu} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"mftMuonProductionVzVsProcess", "MFT-labeled #mu^{#pm} production z by muon process;z_{#mu} (cm);muon TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"mftPiKProductionVz", "MFT-labeled #pi/K production z for #pi/K with direct #mu daughter;z_{#pi/K} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"mftPiKProductionVzVsProcess", "MFT-labeled #pi/K production z by #pi/K process for #pi/K with direct #mu daughter;z_{#pi/K} (cm);#pi/K TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
      {"mftPiKDecayMuonVz", "Direct #mu daughter production z from MFT-labeled #pi/K;z_{#mu} (cm);counts", {HistType::kTH1F, {cfgVzAxis}}},
      {"mftPiKDecayMuonVzVsPiKProcess", "Direct #mu daughter production z by MFT-labeled #pi/K process;z_{#mu} (cm);#pi/K TMCProcess;counts", {HistType::kTH2F, {cfgVzAxis, cfgProcessAxis}}},
    },
    OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    auto counts = registry.get<TH1>(HIST("muonTrackCounts"));
    counts->GetXaxis()->SetBinLabel(1, "generated #mu");
    counts->GetXaxis()->SetBinLabel(2, "MCH-MID #mu tracks");
    counts->GetXaxis()->SetBinLabel(3, "MFT #mu tracks");
    counts->GetXaxis()->SetBinLabel(4, "MCH-MID #mu without MFT");
    counts->GetXaxis()->SetBinLabel(5, "MFT #mu without MCH-MID");
  }

  // Basic MC-particle selections.
  template <typename TMcParticle>
  bool isInConfiguredEtaRange(const TMcParticle& mcParticle) const
  {
    const float etaMin = std::min(cfgEtaMin.value, cfgEtaMax.value);
    const float etaMax = std::max(cfgEtaMin.value, cfgEtaMax.value);
    const float eta = mcParticle.eta();
    return std::isfinite(eta) && eta > etaMin && eta < etaMax;
  }

  template <typename TMcParticle>
  bool isSelectedMuonMcParticle(const TMcParticle& mcParticle) const
  {
    return std::abs(mcParticle.pdgCode()) == 13 && isInConfiguredEtaRange(mcParticle);
  }

  bool isPiK(int pdg) const
  {
    const int absPdg = std::abs(pdg);
    return absPdg == 211 || absPdg == 321 || absPdg == 130 || absPdg == 310;
  }

  template <typename TMcParticle>
  bool hasPiKMother(const TMcParticle& mcParticle) const
  {
    for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
      if (isPiK(mother.pdgCode())) {
        return true;
      }
    }
    return false;
  }

  template <typename TMcParticle>
  void fillMuonVzByMotherProcess(const TMcParticle& mcParticle)
  {
    for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
      registry.fill(HIST("muonVzVsMotherProcess"), mcParticle.vz(), mother.getProcess());
      registry.fill(HIST("motherVzVsMotherProcess"), mother.vz(), mother.getProcess());
    }
  }

  template <typename TMcParticle>
  bool hasFiniteVz(const TMcParticle& mcParticle) const
  {
    return std::isfinite(mcParticle.vz());
  }

  float wrappedDeltaPhi(float phiMuon, float phiMother) const
  {
    const float deltaPhi = phiMuon - phiMother;
    return std::atan2(std::sin(deltaPhi), std::cos(deltaPhi));
  }

  template <typename TMcParticle, typename TMother>
  void fillPiKToMuonDelta(const TMcParticle& mcParticle, const TMother& mother)
  {
    const float deltaPt = mcParticle.pt() - mother.pt();
    const float deltaEta = mcParticle.eta() - mother.eta();
    const float deltaPhi = wrappedDeltaPhi(mcParticle.phi(), mother.phi());
    const float deltaVx = mcParticle.vx() - mother.vx();
    const float deltaVy = mcParticle.vy() - mother.vy();
    const float deltaVz = mcParticle.vz() - mother.vz();
    const float decayLength = std::sqrt(deltaVx * deltaVx + deltaVy * deltaVy + deltaVz * deltaVz);
    if (!std::isfinite(deltaPt) || !std::isfinite(deltaEta) || !std::isfinite(deltaPhi) || !std::isfinite(decayLength)) {
      return;
    }

    registry.fill(HIST("piKToMuonDeltaPt"), deltaPt);
    registry.fill(HIST("piKToMuonDeltaEta"), deltaEta);
    registry.fill(HIST("piKToMuonDeltaPhi"), deltaPhi);
    registry.fill(HIST("piKToMuonDeltaVz"), deltaVz);
    registry.fill(HIST("piKToMuonDecayLength"), decayLength);
    registry.fill(HIST("piKToMuonDecayLengthVsMuonVz"), mcParticle.vz(), decayLength);

    const int motherPdg = std::abs(mother.pdgCode());
    if (motherPdg == 211) {
      registry.fill(HIST("piToMuonDeltaPt"), deltaPt);
      registry.fill(HIST("piToMuonDeltaEta"), deltaEta);
      registry.fill(HIST("piToMuonDeltaPhi"), deltaPhi);
      registry.fill(HIST("piToMuonDeltaVz"), deltaVz);
      registry.fill(HIST("piToMuonDecayLength"), decayLength);
      registry.fill(HIST("piToMuonDecayLengthVsMuonVz"), mcParticle.vz(), decayLength);
    } else {
      registry.fill(HIST("kToMuonDeltaPt"), deltaPt);
      registry.fill(HIST("kToMuonDeltaEta"), deltaEta);
      registry.fill(HIST("kToMuonDeltaPhi"), deltaPhi);
      registry.fill(HIST("kToMuonDeltaVz"), deltaVz);
      registry.fill(HIST("kToMuonDecayLength"), decayLength);
      registry.fill(HIST("kToMuonDecayLengthVsMuonVz"), mcParticle.vz(), decayLength);
    }
  }

  // Processes treated as material-interaction origins. A muon is also tagged as material-origin
  // when any particle in its mother chain has one of these production processes.
  bool isMaterialInteractionProcess(int process) const
  {
    return process == kPHadronic ||
           process == kPHElastic ||
           process == kPHIElastic ||
           process == kPHCElastic ||
           process == kPHInhelastic ||
           process == kPPhotonInhelastic ||
           process == kPEvaporation ||
           process == kPNuclearFission ||
           process == kPNuclearAbsorption ||
           process == kPPbarAnnihilation ||
           process == kPNbarAnnihilation ||
           process == kPNCapture ||
           process == kPMuonNuclear ||
           process == kPElectronNuclear ||
           process == kPPositronNuclear ||
           process == kPPhotoNuclear ||
           process == kPPair;
  }

  template <typename TMcParticle>
  bool hasMaterialInteractionAncestorImpl(const TMcParticle& mcParticle, std::unordered_set<int64_t>& visited) const
  {
    if (!visited.insert(mcParticle.globalIndex()).second) {
      return false;
    }
    for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
      if (isMaterialInteractionProcess(mother.getProcess())) {
        return true;
      }
      if (hasMaterialInteractionAncestorImpl(mother, visited)) {
        return true;
      }
    }
    return false;
  }

  template <typename TMcParticle>
  bool hasMaterialInteractionAncestor(const TMcParticle& mcParticle) const
  {
    std::unordered_set<int64_t> visited;
    return hasMaterialInteractionAncestorImpl(mcParticle, visited);
  }

  template <typename TMcParticle>
  void fillPiKDecayDiagnostics(const TMcParticle& mcParticle)
  {
    for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
      if (!isPiK(mother.pdgCode())) {
        continue;
      }

      const float muonVz = mcParticle.vz();
      const float motherVz = mother.vz();
      fillPiKToMuonDelta(mcParticle, mother);
      if (mother.producedByGenerator()) {
        registry.fill(HIST("forwardMuonVzPiKFromGenerator"), muonVz);
        registry.fill(HIST("piKMotherVzFromGenerator"), motherVz);
        return;
      }

      if (isMaterialInteractionProcess(mother.getProcess()) || hasMaterialInteractionAncestor(mother)) {
        registry.fill(HIST("forwardMuonVzPiKFromMaterial"), muonVz);
        registry.fill(HIST("piKMotherVzFromMaterial"), motherVz);
        return;
      }

      if (mother.getProcess() == kPDecay) {
        registry.fill(HIST("forwardMuonVzPiKFromDecay"), muonVz);
        registry.fill(HIST("piKMotherVzFromDecay"), motherVz);
        return;
      }

      registry.fill(HIST("forwardMuonVzPiKFromOther"), muonVz);
      registry.fill(HIST("piKMotherVzFromOther"), motherVz);
      return;
    }
  }

  // Classify generated forward muons by production mechanism. Pi/K decay diagnostics are filled
  // for all direct pi/K -> mu decays before the muon is assigned to the exclusive Vz categories.
  template <typename TMcParticle>
  void fillGeneratedMuonVzCategories(const TMcParticle& mcParticle)
  {
    const float vz = mcParticle.vz();
    if (mcParticle.producedByGenerator()) {
      registry.fill(HIST("forwardMuonVzGenerator"), vz);
      return;
    }

    const int process = mcParticle.getProcess();
    const bool isPiKDecayMuon = process == kPDecay && hasPiKMother(mcParticle);
    if (isPiKDecayMuon) {
      fillPiKDecayDiagnostics(mcParticle);
    }

    if (isMaterialInteractionProcess(process) || hasMaterialInteractionAncestor(mcParticle)) {
      registry.fill(HIST("forwardMuonVzMaterial"), vz);
      return;
    }

    if (process == kPDecay) {
      registry.fill(HIST("forwardMuonVzNaturalDecay"), vz);
      if (hasPiKMother(mcParticle)) {
        registry.fill(HIST("forwardMuonVzPiKDecay"), vz);
      }
      return;
    }

    registry.fill(HIST("forwardMuonVzOther"), vz);
  }

  template <typename TMcParticles>
  void fillGeneratedMuonObservables(const TMcParticles& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (!isSelectedMuonMcParticle(mcParticle)) {
        continue;
      }
      registry.fill(HIST("forwardMuonVz"), mcParticle.vz());
      registry.fill(HIST("muonVzVsMuonProcess"), mcParticle.vz(), mcParticle.getProcess());
      fillMuonVzByMotherProcess(mcParticle);
      fillGeneratedMuonVzCategories(mcParticle);
      registry.fill(HIST("muonTrackCounts"), 1.);
    }
  }

  template <typename TMcParticles>
  void fillPiKMotherDaughterObservables(const TMcParticles& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (!std::isfinite(mcParticle.vz())) {
        continue;
      }

      for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
        const int absMotherPdg = std::abs(mother.pdgCode());
        if (!isPiK(absMotherPdg)) {
          continue;
        }

        const int process = mcParticle.getProcess();
        const int pdg = mcParticle.pdgCode();
        const bool isMuon = std::abs(pdg) == 13;
        registry.fill(HIST("piKMotherDaughterVz"), mcParticle.vz());
        registry.fill(HIST("piKMotherDaughterVzVsProcess"), mcParticle.vz(), process);
        if (!isMuon) {
          registry.fill(HIST("piKMotherNonMuonDaughterVz"), mcParticle.vz());
          registry.fill(HIST("piKMotherNonMuonDaughterVzVsProcess"), mcParticle.vz(), process);
          registry.fill(HIST("piKMotherNonMuonDaughterVzVsPdg"), mcParticle.vz(), pdg);
        }
        if (absMotherPdg == 211) {
          registry.fill(HIST("piMotherDaughterVz"), mcParticle.vz());
          registry.fill(HIST("piMotherDaughterVzVsProcess"), mcParticle.vz(), process);
          if (!isMuon) {
            registry.fill(HIST("piMotherNonMuonDaughterVz"), mcParticle.vz());
            registry.fill(HIST("piMotherNonMuonDaughterVzVsProcess"), mcParticle.vz(), process);
            registry.fill(HIST("piMotherNonMuonDaughterVzVsPdg"), mcParticle.vz(), pdg);
          }
        } else {
          registry.fill(HIST("kMotherDaughterVz"), mcParticle.vz());
          registry.fill(HIST("kMotherDaughterVzVsProcess"), mcParticle.vz(), process);
          if (!isMuon) {
            registry.fill(HIST("kMotherNonMuonDaughterVz"), mcParticle.vz());
            registry.fill(HIST("kMotherNonMuonDaughterVzVsProcess"), mcParticle.vz(), process);
            registry.fill(HIST("kMotherNonMuonDaughterVzVsPdg"), mcParticle.vz(), pdg);
          }
        }
      }
    }
  }

  template <typename TFwdTracks>
  void collectMchMidMuonMcIds(const TFwdTracks& fwdTracks, std::unordered_set<int64_t>& selectedMchMidMuonMcIds)
  {
    for (const auto& fwdTrack : fwdTracks) {
      if (fwdTrack.trackType() != 3 || !fwdTrack.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = fwdTrack.mcParticle();
      if (!isSelectedMuonMcParticle(mcParticle)) {
        continue;
      }

      selectedMchMidMuonMcIds.insert(mcParticle.globalIndex());
      registry.fill(HIST("muonTrackCounts"), 2.);
    }
  }

  template <typename TFwdTracks>
  void fillFwdTrackMuonPurity(const TFwdTracks& fwdTracks)
  {
    for (const auto& fwdTrack : fwdTracks) {
      if (!fwdTrack.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = fwdTrack.mcParticle();
      if (!isInConfiguredEtaRange(mcParticle) || !std::isfinite(mcParticle.pt())) {
        continue;
      }

      const bool isMuon = std::abs(mcParticle.pdgCode()) == 13;
      registry.fill(HIST("fwdTrackTruePtAll"), mcParticle.pt());
      registry.fill(HIST("fwdTrackMuonFractionVsTruePt"), mcParticle.pt(), isMuon ? 1. : 0.);
      if (isMuon) {
        registry.fill(HIST("fwdTrackTruePtMuon"), mcParticle.pt());
      }
    }
  }

  template <typename TMftTracks>
  void collectMftMuonMcIds(const TMftTracks& mftTracks, std::unordered_set<int64_t>& selectedMftMuonMcIds)
  {
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = mftTrack.mcParticle();
      if (!isSelectedMuonMcParticle(mcParticle)) {
        continue;
      }

      selectedMftMuonMcIds.insert(mcParticle.globalIndex());
      registry.fill(HIST("muonTrackCounts"), 3.);
    }
  }

  template <typename TMcParticle>
  void fillMftMuonLabelProductionVz(const TMcParticle& mcParticle)
  {
    if (!hasFiniteVz(mcParticle)) {
      return;
    }
    registry.fill(HIST("mftMuonProductionVz"), mcParticle.vz());
    registry.fill(HIST("mftMuonProductionVzVsProcess"), mcParticle.vz(), mcParticle.getProcess());
  }

  template <typename TMcParticle>
  void fillMftPiKLabelProductionAndDecayVz(const TMcParticle& mcParticle)
  {
    if (!mcParticle.has_daughters()) {
      return;
    }

    bool hasDirectMuonDaughter = false;
    for (const auto& daughter : mcParticle.template daughters_as<aod::McParticles>()) {
      if (std::abs(daughter.pdgCode()) != 13 || !hasFiniteVz(daughter)) {
        continue;
      }

      hasDirectMuonDaughter = true;
      registry.fill(HIST("mftPiKDecayMuonVz"), daughter.vz());
      registry.fill(HIST("mftPiKDecayMuonVzVsPiKProcess"), daughter.vz(), mcParticle.getProcess());
    }

    if (hasDirectMuonDaughter && hasFiniteVz(mcParticle)) {
      registry.fill(HIST("mftPiKProductionVz"), mcParticle.vz());
      registry.fill(HIST("mftPiKProductionVzVsProcess"), mcParticle.vz(), mcParticle.getProcess());
    }
  }

  template <typename TMftTracks>
  void fillMftMcLabelVzDiagnostics(const TMftTracks& mftTracks)
  {
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_mcParticle()) {
        continue;
      }

      const auto& mcParticle = mftTrack.mcParticle();
      const int absPdg = std::abs(mcParticle.pdgCode());
      if (absPdg == 13) {
        fillMftMuonLabelProductionVz(mcParticle);
      } else if (isPiK(absPdg)) {
        fillMftPiKLabelProductionAndDecayVz(mcParticle);
      }
    }
  }

  template <typename TFwdTracks>
  void fillMchMidMuonsWithoutMft(const TFwdTracks& fwdTracks, const std::unordered_set<int64_t>& selectedMftMuonMcIds)
  {
    for (const auto& fwdTrack : fwdTracks) {
      if (fwdTrack.trackType() != 3 || !fwdTrack.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = fwdTrack.mcParticle();
      if (!isSelectedMuonMcParticle(mcParticle)) {
        continue;
      }
      if (selectedMftMuonMcIds.find(mcParticle.globalIndex()) == selectedMftMuonMcIds.end()) {
        registry.fill(HIST("muonTrackCounts"), 4.);
      }
    }
  }

  template <typename TMftTracks>
  void fillMftMuonsWithoutMchMid(const TMftTracks& mftTracks, const std::unordered_set<int64_t>& selectedMchMidMuonMcIds)
  {
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = mftTrack.mcParticle();
      if (!isSelectedMuonMcParticle(mcParticle)) {
        continue;
      }
      if (selectedMchMidMuonMcIds.find(mcParticle.globalIndex()) == selectedMchMidMuonMcIds.end()) {
        registry.fill(HIST("muonTrackCounts"), 5.);
      }
    }
  }

  void process(aod::McParticles const& mcParticles,
               LabeledFwdTracks const& fwdTracks,
               LabeledMftTracks const& mftTracks)
  {
    std::unordered_set<int64_t> selectedMchMidMuonMcIds;
    std::unordered_set<int64_t> selectedMftMuonMcIds;

    fillGeneratedMuonObservables(mcParticles);
    fillPiKMotherDaughterObservables(mcParticles);
    collectMchMidMuonMcIds(fwdTracks, selectedMchMidMuonMcIds);
    fillFwdTrackMuonPurity(fwdTracks);
    collectMftMuonMcIds(mftTracks, selectedMftMuonMcIds);
    fillMftMcLabelVzDiagnostics(mftTracks);
    fillMchMidMuonsWithoutMft(fwdTracks, selectedMftMuonMcIds);
    fillMftMuonsWithoutMchMid(mftTracks, selectedMchMidMuonMcIds);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<McForwardMuonVz>(cfgc, TaskName{"mc-forward-muon-vz"})};
}
