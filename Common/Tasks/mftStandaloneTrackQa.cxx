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

#include "Common/DataModel/CollisionAssociationTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TAxis.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct MFTStandaloneTrackQa {
  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>;
  using CollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using MFTTracksMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<std::string> grpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  int runNumber = -1;
  float bzAtMFTCenter = 0.f;
  static constexpr double centerMFT[3] = {0., 0., -61.4};

  struct : ConfigurableGroup {
    Configurable<float> etaMin{"cfgEtaMin", -3.6f, "Minimum track/truth eta (exclusive)"};
    Configurable<float> etaMax{"cfgEtaMax", -2.5f, "Maximum track/truth eta (exclusive)"};
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
  ConfigurableAxis axisDcaX{"axisDcaX", {400, -20., 20.}, "#it{DCA}_{x} (cm)"};
  ConfigurableAxis axisDcaY{"axisDcaY", {400, -20., 20.}, "#it{DCA}_{y} (cm)"};
  ConfigurableAxis axisDcaXY{"axisDcaXY", {400, 0., 20.}, "#it{DCA}_{xy} (cm)"};
  ConfigurableAxis axisBestDcaXY{"axisBestDcaXY", {400, 0., 20.}, "best #it{DCA}_{xy} (cm)"};
  ConfigurableAxis axisNonBestCollisionDeltaTime{"axisNonBestCollisionDeltaTime", {2000, -25000., 25000.}, "#Delta #it{t}_{other-best} (ns)"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {5000, -0.5, 4999.5}, "selected MFT tracks per collision"};
  ConfigurableAxis axisMultiplicityCoarse{"axisMultiplicityCoarse", {200, -0.5, 4999.5}, "selected MFT tracks per collision"};
  ConfigurableAxis axisPtRelativeResidual{"axisPtRelativeResidual", {400, -5., 5.}, "(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true}"};
  ConfigurableAxis axisEtaResidual{"axisEtaResidual", {400, -0.2, 0.2}, "#it{#eta}^{reco} - #it{#eta}^{true}"};
  ConfigurableAxis axisPhiResidual{"axisPhiResidual", {400, -0.2, 0.2}, "#it{#varphi}^{reco} - #it{#varphi}^{true} (rad)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Category prefixes are part of every histogram name because the analysis
  // output sink merges objects by task and histogram name, not by registry name.
  HistogramRegistry registryAssociatedUnique{"associatedUnique", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryAssociatedAmbiguous{"associatedAmbiguous", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryAssociatedAmbiguousBest{"associatedAmbiguousBest", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryAssociatedUnassigned{"associatedUnassigned", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryAssociationSummary{"association", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry* activeRegistry = nullptr;
  std::string histogramBookingPrefix;

  struct TrackCounts {
    int all = 0;
    int ltf = 0;
    int ca = 0;
  };

  struct PropagatedTrackValues {
    float x = 0.f;
    float y = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    float pt = 0.f;
    float dcaX = 0.f;
    float dcaY = 0.f;
    float dcaXY = 0.f;
  };

  struct BestCollisionResult {
    bool found = false;
    int64_t collisionId = -1;
    float dcaXY = std::numeric_limits<float>::max();
    PropagatedTrackValues propagatedValues;
  };

  struct AssociationTruthSummary {
    bool hasTruthMatchedCandidate = false;
    bool hasCandidateDcaXY = false;
    float candidateDcaXY = 0.f;
  };

  enum RecoMask : uint8_t {
    kLTF = 0x1,
    kCA = 0x2
  };

  enum class AssociationCategory : uint8_t {
    Unassigned = 0,
    Unique = 1,
    Ambiguous = 2
  };

  enum class AssociationTruthStatus : uint8_t {
    TrackTruthUnavailable = 0,
    TruthCollisionNotReconstructed = 1,
    TruthCollisionReconstructedAssociationCorrect = 2,
    TruthCollisionReconstructedAssociationWrong = 3,
    TruthCollisionSplitAssociationCorrect = 4,
    TruthCollisionSplitAssociationWrong = 5
  };

  using AssociationCounts = std::unordered_map<int64_t, int>;
  using BestCollisionResults = std::unordered_map<int64_t, BestCollisionResult>;
  using AssociationTruthSummaries = std::unordered_map<int64_t, AssociationTruthSummary>;
  using ReconstructedCollisionCounts = std::unordered_map<int64_t, int>;

  std::string prefixedHistogramName(char const* relativeName) const
  {
    return histogramBookingPrefix + relativeName;
  }

  HistPtr addHistogram(char const* name, char const* title, HistType type, std::vector<AxisSpec> const& axes, bool callSumw2 = false)
  {
    return activeRegistry->add(prefixedHistogramName(name), title, type, axes, callSumw2);
  }

  template <typename T>
  std::shared_ptr<T> addHistogram(char const* name, char const* title, HistType type, std::vector<AxisSpec> const& axes, bool callSumw2 = false)
  {
    return activeRegistry->add<T>(prefixedHistogramName(name), title, type, axes, callSumw2);
  }

  void addHistogramClone(char const* source, char const* target)
  {
    activeRegistry->addClone(prefixedHistogramName(source), prefixedHistogramName(target));
  }

  template <typename TName, typename... Ts>
  void fillHistogram(TName const& name, Ts... values)
  {
    if (activeRegistry == &registryAssociatedUnique) {
      activeRegistry->fill(HIST("associatedUnique/") + name, values...);
    } else if (activeRegistry == &registryAssociatedAmbiguous) {
      activeRegistry->fill(HIST("associatedAmbiguous/") + name, values...);
    } else if (activeRegistry == &registryAssociatedAmbiguousBest) {
      activeRegistry->fill(HIST("associatedAmbiguousBest/") + name, values...);
    } else if (activeRegistry == &registryAssociatedUnassigned) {
      activeRegistry->fill(HIST("associatedUnassigned/") + name, values...);
    } else {
      activeRegistry->fill(name, values...);
    }
  }

  template <typename T, typename TName>
  std::shared_ptr<T> getHistogram(TName const& name)
  {
    if (activeRegistry == &registryAssociatedUnique) {
      return activeRegistry->get<T>(HIST("associatedUnique/") + name);
    }
    if (activeRegistry == &registryAssociatedAmbiguous) {
      return activeRegistry->get<T>(HIST("associatedAmbiguous/") + name);
    }
    if (activeRegistry == &registryAssociatedAmbiguousBest) {
      return activeRegistry->get<T>(HIST("associatedAmbiguousBest/") + name);
    }
    if (activeRegistry == &registryAssociatedUnassigned) {
      return activeRegistry->get<T>(HIST("associatedUnassigned/") + name);
    }
    return activeRegistry->get<T>(name);
  }

  void addRecoHistograms()
  {
    addHistogram("reco/all/pt", "MFT standalone tracks;#it{p}_{T} (GeV/#it{c});counts", HistType::kTH1F, {axisPt});
    addHistogram("reco/all/eta", "MFT standalone tracks;#it{#eta};counts", HistType::kTH1F, {axisEta});
    addHistogram("reco/all/phi", "MFT standalone tracks;#it{#varphi} (rad);counts", HistType::kTH1F, {axisPhi});
    addHistogram("reco/all/x", "MFT standalone tracks at stored track #it{z};#it{x} (cm);counts", HistType::kTH1F, {axisXY});
    addHistogram("reco/all/y", "MFT standalone tracks at stored track #it{z};#it{y} (cm);counts", HistType::kTH1F, {axisXY});
    addHistogram("reco/all/z", "MFT standalone tracks;stored track #it{z} (cm);counts", HistType::kTH1F, {axisZ});
    addHistogram("reco/all/signed1Pt", "MFT standalone tracks;#it{q}/#it{p}_{T} (#it{c}/GeV);counts", HistType::kTH1F, {axisSigned1Pt});
    addHistogram("reco/all/sign", "MFT standalone tracks;charge sign;counts", HistType::kTH1F, {{3, -1.5, 1.5}});
    addHistogram("reco/all/nClusters", "MFT standalone tracks;MFT clusters;counts", HistType::kTH1F, {axisNClusters});
    addHistogram("reco/all/chi2Ndf", "MFT standalone tracks;#chi^{2}/NDF;counts", HistType::kTH1F, {axisChi2Ndf});
    addHistogram("reco/all/xy", "MFT standalone tracks at stored track #it{z};#it{x} (cm);#it{y} (cm)", HistType::kTH2F, {axisXY, axisXY});
    addHistogram("reco/all/etaPhi", "MFT standalone tracks;#it{#eta};#it{#varphi} (rad)", HistType::kTH2F, {axisEta, axisPhi});
    addHistogram("reco/all/nClustersVsEta", "MFT standalone tracks;#it{#eta};MFT clusters", HistType::kTH2F, {axisEta, axisNClusters});
    addHistogram("reco/all/chi2NdfVsPt", "MFT standalone tracks;#it{p}_{T} (GeV/#it{c});#chi^{2}/NDF", HistType::kTH2F, {axisPt, axisChi2Ndf});
    addHistogram("reco/all/trackTime", "MFT standalone tracks;track time (ns);counts", HistType::kTH1F, {axisTrackTime});
    addHistogram("reco/all/trackTimeRes", "MFT standalone tracks;track-time resolution (ns);counts", HistType::kTH1F, {axisTrackTimeRes});
    addHistogram("reco/all/layerOccupancy", "MFT clusters used by standalone tracks;MFT layer;clusters", HistType::kTH1F, {{10, -0.5, 9.5}});
    addHistogram("reco/all/clusterSizeVsLayer", "MFT cluster size by layer;MFT layer;cluster size", HistType::kTH2F, {{10, -0.5, 9.5}, axisClusterSize});
    addHistogram("reco/all/layerPattern", "MFT standalone-track layer pattern;10-layer pattern;tracks", HistType::kTH1F, {{1024, -0.5, 1023.5}});
    addHistogramClone("reco/all/", "reco/ltf/");
    addHistogramClone("reco/all/", "reco/ca/");

    addHistogramClone("reco/all/", "reco/selected/all/");
    addHistogramClone("reco/selected/all/", "reco/selected/ltf/");
    addHistogramClone("reco/selected/all/", "reco/selected/ca/");
  }

  void addPropagatedHistograms()
  {
    addHistogram("propagated/pv/status", "MFT track helix propagation to associated PV #it{z} status;;tracks", HistType::kTH1F, {{3, -0.5, 2.5}});
    addHistogram("propagated/pv/all/pt", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{p}_{T} (GeV/#it{c});counts", HistType::kTH1F, {axisPt});
    addHistogram("propagated/pv/all/eta", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{#eta};counts", HistType::kTH1F, {axisEta});
    addHistogram("propagated/pv/all/phi", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{#varphi} (rad);counts", HistType::kTH1F, {axisPhi});
    addHistogram("propagated/pv/all/x", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{x} (cm);counts", HistType::kTH1F, {axisXY});
    addHistogram("propagated/pv/all/y", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{y} (cm);counts", HistType::kTH1F, {axisXY});
    addHistogram("propagated/pv/all/xy", "MFT standalone tracks helix-propagated to associated PV #it{z};#it{x} (cm);#it{y} (cm)", HistType::kTH2F, {axisXY, axisXY});
    addHistogram("propagated/pv/all/dcaX", "MFT standalone tracks helix-propagated to associated PV;#it{DCA}_{x} (cm);counts", HistType::kTH1F, {axisDcaX});
    addHistogram("propagated/pv/all/dcaY", "MFT standalone tracks helix-propagated to associated PV;#it{DCA}_{y} (cm);counts", HistType::kTH1F, {axisDcaY});
    addHistogram("propagated/pv/all/dcaXY", "MFT standalone tracks helix-propagated to associated PV;#it{DCA}_{xy} (cm);counts", HistType::kTH1F, {axisDcaXY});
    addHistogramClone("propagated/pv/all/", "propagated/pv/ltf/");
    addHistogramClone("propagated/pv/all/", "propagated/pv/ca/");
  }

  void addCAContributionHistograms()
  {
    addHistogram("caContribution/fractionVsPt", "CA contribution among selected tracks;#it{p}_{T}^{reco} (GeV/#it{c});CA fraction", HistType::kTProfile, {axisPt});
    addHistogram("caContribution/fractionVsEta", "CA contribution among selected tracks;#it{#eta};CA fraction", HistType::kTProfile, {axisEta});
    addHistogram("caContribution/fractionVsNClusters", "CA contribution among selected tracks;MFT clusters;CA fraction", HistType::kTProfile, {axisNClusters});
    addHistogram("caContribution/fractionVsChi2Ndf", "CA contribution among selected tracks;#chi^{2}/NDF;CA fraction", HistType::kTProfile, {axisChi2Ndf});

    addHistogram("collision/nTracksAll", "Selected MFT tracks per collision;MFT tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    addHistogram("collision/nTracksLTF", "Selected LTF tracks per collision;LTF tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    addHistogram("collision/nTracksCA", "Selected CA tracks per collision;CA tracks;collisions", HistType::kTH1F, {axisMultiplicity});
    addHistogram("collision/nCAVsNLTF", "CA and LTF multiplicities;LTF tracks;CA tracks", HistType::kTH2F, {axisMultiplicityCoarse, axisMultiplicityCoarse});
    addHistogram("collision/caFraction", "CA contribution per collision;CA/(LTF+CA);collisions", HistType::kTH1F, {{110, -0.05, 1.05}});
    addHistogram("collision/caFractionVsNTracks", "CA contribution vs selected multiplicity;selected MFT tracks;CA/(LTF+CA)", HistType::kTProfile, {axisMultiplicity});
  }

  void addMCHistograms()
  {
    addHistogram("mc/labelStatusVsAlgorithm", "MC label status by tracking algorithm;tracking algorithm;MC label status", HistType::kTH2F, {{2, -0.5, 1.5}, {5, -0.5, 4.5}});

    addHistogram("mc/ltf/nonCorrectLabelFractionVsPt", "Selected MFT tracks;#it{p}_{T}^{reco} (GeV/#it{c});non-correct label fraction", HistType::kTProfile, {axisPt});
    addHistogram("mc/ltf/nonCorrectLabelFractionVsEta", "Selected MFT tracks;#it{#eta}^{reco};non-correct label fraction", HistType::kTProfile, {axisEta});
    addHistogram("mc/ltf/nonCorrectLabelFractionVsNClusters", "Selected MFT tracks;MFT clusters;non-correct label fraction", HistType::kTProfile, {axisNClusters});
    addHistogram("mc/ltf/ptRelativeResidualVsTruePt", "Selected truth-matched MFT tracks;#it{p}_{T}^{true} (GeV/#it{c});(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true}", HistType::kTH2F, {axisPt, axisPtRelativeResidual});
    addHistogram("mc/ltf/etaResidualVsTrueEta", "Selected truth-matched MFT tracks;#it{#eta}^{true};#it{#eta}^{reco} - #it{#eta}^{true}", HistType::kTH2F, {axisEta, axisEtaResidual});
    addHistogram("mc/ltf/phiResidualVsTrueEta", "Selected truth-matched MFT tracks;#it{#eta}^{true};#it{#varphi}^{reco} - #it{#varphi}^{true} (rad)", HistType::kTH2F, {axisEta, axisPhiResidual});
    addHistogram("mc/ltf/chargeCorrectVsTruePt", "Selected truth-matched MFT tracks;#it{p}_{T}^{true} (GeV/#it{c});correct charge fraction", HistType::kTProfile, {axisPt});
    addHistogramClone("mc/ltf/", "mc/ca/");

    addHistogram("mc/truth/denominator/pt", "Selected charged truth particles;#it{p}_{T}^{true} (GeV/#it{c});particles", HistType::kTH1F, {axisPt});
    addHistogram("mc/truth/denominator/eta", "Selected charged truth particles;#it{#eta}^{true};particles", HistType::kTH1F, {axisEta});
    addHistogram("mc/truth/denominator/ptEta", "Selected charged truth particles;#it{p}_{T}^{true} (GeV/#it{c});#it{#eta}^{true}", HistType::kTH2F, {axisPt, axisEta});
    addHistogramClone("mc/truth/denominator/", "mc/truth/inclusive/");
    addHistogramClone("mc/truth/denominator/", "mc/truth/ltf/");
    addHistogramClone("mc/truth/denominator/", "mc/truth/ca/");
    addHistogramClone("mc/truth/denominator/", "mc/truth/ltfOnly/");
    addHistogramClone("mc/truth/denominator/", "mc/truth/caOnly/");
    addHistogramClone("mc/truth/denominator/", "mc/truth/both/");
    addHistogram("mc/truth/recoveryCategory", "Selected truth-particle recovery category;recovery category;particles", HistType::kTH1F, {{4, -0.5, 3.5}});
    addHistogram("mc/truth/nRecoTracks", "Selected reconstructed tracks per truth particle;reconstructed tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
    addHistogram("mc/truth/nRecoTracksLTF", "Selected LTF tracks per truth particle;LTF tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
    addHistogram("mc/truth/nRecoTracksCA", "Selected CA tracks per truth particle;CA tracks;particles", HistType::kTH1F, {{11, -0.5, 10.5}});
  }

  void setAssociationTruthStatusLabels(TAxis* axis) const
  {
    axis->SetBinLabel(1, "track truth unavailable");
    axis->SetBinLabel(2, "truth collision not reconstructed");
    axis->SetBinLabel(3, "single reco collision: association correct");
    axis->SetBinLabel(4, "single reco collision: association wrong");
    axis->SetBinLabel(5, "split collision: association correct");
    axis->SetBinLabel(6, "split collision: association wrong");
  }

  void addAssociationTruthHistograms()
  {
    auto statusAll = addHistogram<TH1>("mc/association/truthStatus/all", "Hierarchical MC collision-reconstruction and MFT-association status;truth status;MFT tracks", HistType::kTH1F, {{6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusAll->GetXaxis());
    auto statusSelected = addHistogram<TH1>("mc/association/truthStatus/selected", "Hierarchical MC collision-reconstruction and MFT-association status for selected tracks;truth status;MFT tracks", HistType::kTH1F, {{6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusSelected->GetXaxis());
    auto statusVsPt = addHistogram<TH2>("mc/association/truthStatusVsPt", "Hierarchical MC collision-reconstruction and MFT-association status for selected tracks;#it{p}_{T}^{reco} (GeV/#it{c});truth status", HistType::kTH2F, {axisPt, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusVsPt->GetYaxis());
    auto statusVsEta = addHistogram<TH2>("mc/association/truthStatusVsEta", "Hierarchical MC collision-reconstruction and MFT-association status for selected tracks;#it{#eta}^{reco};truth status", HistType::kTH2F, {axisEta, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusVsEta->GetYaxis());
    auto statusVsPhi = addHistogram<TH2>("mc/association/truthStatusVsPhi", "Hierarchical MC collision-reconstruction and MFT-association status for selected tracks;#it{#varphi}^{reco} (rad);truth status", HistType::kTH2F, {axisPhi, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusVsPhi->GetYaxis());
    auto statusVsChi2Ndf = addHistogram<TH2>("mc/association/truthStatusVsChi2Ndf", "Hierarchical MC collision-reconstruction and MFT-association status for selected tracks;MFT tracking #chi^{2}/NDF;truth status", HistType::kTH2F, {axisChi2Ndf, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(statusVsChi2Ndf->GetYaxis());
  }

  void initializeAssociationTruthHistograms()
  {
    activeRegistry = &registryAssociatedUnique;
    histogramBookingPrefix = "associatedUnique/";
    addAssociationTruthHistograms();
    auto uniqueStatusVsDcaXY = addHistogram<TH2>("mc/association/truthStatusVsDcaXY", "Hierarchical MC truth status of selected uniquely associated tracks with successful propagation;#it{DCA}_{xy} (cm);truth status", HistType::kTH2F, {axisDcaXY, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(uniqueStatusVsDcaXY->GetYaxis());

    activeRegistry = &registryAssociatedAmbiguous;
    histogramBookingPrefix = "associatedAmbiguous/";
    addAssociationTruthHistograms();

    activeRegistry = &registryAssociatedAmbiguousBest;
    histogramBookingPrefix = "associatedAmbiguousBest/";
    addAssociationTruthHistograms();
    auto bestStatusVsDcaXY = addHistogram<TH2>("mc/association/truthStatusVsDcaXY", "Hierarchical MC truth status of selected DCA-best associations;best #it{DCA}_{xy} (cm);truth status", HistType::kTH2F, {axisBestDcaXY, {6, -0.5, 5.5}});
    setAssociationTruthStatusLabels(bestStatusVsDcaXY->GetYaxis());

    activeRegistry = &registryAssociatedUnassigned;
    histogramBookingPrefix = "associatedUnassigned/";
    addAssociationTruthHistograms();
  }

  void initializeHistogramRegistry(HistogramRegistry& outputRegistry, char const* prefix, bool withMC)
  {
    activeRegistry = &outputRegistry;
    histogramBookingPrefix = prefix;

    addRecoHistograms();
    addPropagatedHistograms();
    addCAContributionHistograms();
    if (withMC) {
      addMCHistograms();
    }

    auto selectionCounter = addHistogram<TH1>("reco/selectionCounter", "MFT standalone track selection;;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    selectionCounter->GetXaxis()->SetBinLabel(1, "all");
    selectionCounter->GetXaxis()->SetBinLabel(2, "selected");

    auto algorithm = addHistogram<TH1>("reco/selectedAlgorithm", "Selected MFT standalone tracks;tracking algorithm;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    algorithm->GetXaxis()->SetBinLabel(1, "LTF");
    algorithm->GetXaxis()->SetBinLabel(2, "CA");

    auto association = addHistogram<TH1>("reco/associationStatus", "MFT standalone track collision association;association;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    association->GetXaxis()->SetBinLabel(1, "unassociated");
    association->GetXaxis()->SetBinLabel(2, "associated");

    auto propagationStatus = getHistogram<TH1>(HIST("propagated/pv/status"));
    propagationStatus->GetXaxis()->SetBinLabel(1, "no collision");
    propagationStatus->GetXaxis()->SetBinLabel(2, "failed");
    propagationStatus->GetXaxis()->SetBinLabel(3, "success");

    if (withMC) {
      auto labelStatus = getHistogram<TH2>(HIST("mc/labelStatusVsAlgorithm"));
      labelStatus->GetXaxis()->SetBinLabel(1, "LTF");
      labelStatus->GetXaxis()->SetBinLabel(2, "CA");
      labelStatus->GetYaxis()->SetBinLabel(1, "correct");
      labelStatus->GetYaxis()->SetBinLabel(2, "fake");
      labelStatus->GetYaxis()->SetBinLabel(3, "noise");
      labelStatus->GetYaxis()->SetBinLabel(4, "fake+noise");
      labelStatus->GetYaxis()->SetBinLabel(5, "no label");

      auto recovery = getHistogram<TH1>(HIST("mc/truth/recoveryCategory"));
      recovery->GetXaxis()->SetBinLabel(1, "missed");
      recovery->GetXaxis()->SetBinLabel(2, "LTF only");
      recovery->GetXaxis()->SetBinLabel(3, "CA only");
      recovery->GetXaxis()->SetBinLabel(4, "LTF and CA");
    }
  }

  void addAmbiguousBestCollisionHistograms()
  {
    auto bestStatus = addHistogram<TH1>("association/bestCollisionStatus", "Best-collision selection for ambiguous MFT tracks;;tracks", HistType::kTH1F, {{2, -0.5, 1.5}});
    bestStatus->GetXaxis()->SetBinLabel(1, "no valid candidate");
    bestStatus->GetXaxis()->SetBinLabel(2, "best found");

    auto candidateStatus = addHistogram<TH1>("association/candidatePropagationStatus", "Propagation of ambiguous MFT tracks to candidate collisions;;track-collision candidates", HistType::kTH1F, {{2, -0.5, 1.5}});
    candidateStatus->GetXaxis()->SetBinLabel(1, "failed");
    candidateStatus->GetXaxis()->SetBinLabel(2, "success");

    addHistogram("association/bestDcaXY", "Best collision selected by minimum transverse distance;best #it{DCA}_{xy} (cm);tracks", HistType::kTH1F, {axisBestDcaXY});
    addHistogram("association/nonBestCollisionDeltaTime/all", "Non-best collision time relative to the best collision;#Delta #it{t}_{other-best} (ns);track-collision candidates", HistType::kTH1F, {axisNonBestCollisionDeltaTime});
    addHistogram("association/nonBestCollisionDeltaTime/selected", "Non-best collision time relative to the best collision for selected tracks;#Delta #it{t}_{other-best} (ns);track-collision candidates", HistType::kTH1F, {axisNonBestCollisionDeltaTime});
  }

  void init(InitContext const&)
  {
    const int nEnabledProcesses = static_cast<int>(doprocessData) + static_cast<int>(doprocessMC) +
                                  static_cast<int>(doprocessDataWithAssociations) + static_cast<int>(doprocessMCWithAssociations);
    if (nEnabledProcesses != 1) {
      LOGF(fatal, "Enable exactly one data/MC process, with or without MFTTrackAssoc");
    }
    if (selection.etaMin.value > selection.etaMax.value || selection.ptMin.value > selection.ptMax.value) {
      LOGF(fatal, "Invalid MFT standalone QA selection range");
    }
    if (selection.minNClusters.value < 0 || selection.maxChi2Ndf.value < 0.f) {
      LOGF(fatal, "MFT cluster and chi2/NDF selections must be non-negative");
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const bool withMC = doprocessMC || doprocessMCWithAssociations;
    const bool withAssociations = doprocessDataWithAssociations || doprocessMCWithAssociations;
    initializeHistogramRegistry(registry, "", withMC);
    if (withAssociations) {
      initializeHistogramRegistry(registryAssociatedUnique, "associatedUnique/", withMC);
      initializeHistogramRegistry(registryAssociatedAmbiguous, "associatedAmbiguous/", withMC);
      initializeHistogramRegistry(registryAssociatedUnassigned, "associatedUnassigned/", withMC);
      initializeHistogramRegistry(registryAssociatedAmbiguousBest, "associatedAmbiguousBest/", withMC);
      addAmbiguousBestCollisionHistograms();
      if (withMC) {
        initializeAssociationTruthHistograms();
      }

      registryAssociationSummary.add("nAssociationsPerTrack", "Collision associations per original MFT track;associations;MFT tracks", HistType::kTH1F, {{21, -0.5, 20.5}});
      auto category = registryAssociationSummary.add<TH1>("trackCategory", "MFT tracks by collision-association category;category;MFT tracks", HistType::kTH1F, {{3, -0.5, 2.5}});
      category->GetXaxis()->SetBinLabel(1, "unassigned");
      category->GetXaxis()->SetBinLabel(2, "unique");
      category->GetXaxis()->SetBinLabel(3, "ambiguous");
    }
    activeRegistry = &registry;
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

  void initMagneticField(ExtBCs const& bcs)
  {
    if (bcs.size() == 0) {
      LOGF(fatal, "Cannot initialize MFT propagation field: BC table is empty");
    }
    const auto bc = bcs.begin();
    if (runNumber == bc.runNumber()) {
      return;
    }

    auto* grpMag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpMagPath, bc.timestamp());
    if (grpMag == nullptr) {
      LOGF(fatal, "Cannot load GRPMagField from %s for run %d", grpMagPath.value.c_str(), bc.runNumber());
    }
    o2::base::Propagator::initFieldFromGRP(grpMag);
    auto* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (field == nullptr) {
      LOGF(fatal, "Magnetic field initialization failed for run %d", bc.runNumber());
    }
    bzAtMFTCenter = field->getBz(centerMFT);
    runNumber = bc.runNumber();
    LOGF(info, "MFT-to-PV propagation uses Bz = %.6f kG sampled at z = %.1f cm for run %d", bzAtMFTCenter, centerMFT[2], runNumber);
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
  bool propagateToPV(TTrack const& track, float pvZ, PropagatedTrackValues& values) const
  {
    if (!std::isfinite(track.x()) || !std::isfinite(track.y()) || !std::isfinite(track.z()) ||
        !std::isfinite(track.phi()) || !std::isfinite(track.tgl()) || !std::isfinite(track.signed1Pt()) ||
        std::abs(track.tgl()) < o2::constants::math::Almost0 ||
        std::abs(track.signed1Pt()) < o2::constants::math::Almost0 || !std::isfinite(pvZ) || !std::isfinite(bzAtMFTCenter)) {
      return false;
    }

    const o2::track::SMatrix5 parameters(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    const o2::track::SMatrix55Sym covariances{};
    // MFT AODs do not store the outward Kalman-refit state. Propagate the stored
    // parameters as a helix without cluster refitting or material corrections.
    o2::track::TrackParCovFwd propagatedTrack(track.z(), parameters, covariances, track.chi2());
    propagatedTrack.propagateToZ(pvZ, bzAtMFTCenter);

    values.x = static_cast<float>(propagatedTrack.getX());
    values.y = static_cast<float>(propagatedTrack.getY());
    values.eta = static_cast<float>(propagatedTrack.getEta());
    values.phi = wrapPhi(static_cast<float>(propagatedTrack.getPhi()));
    values.pt = static_cast<float>(propagatedTrack.getPt());
    return std::isfinite(values.x) && std::isfinite(values.y) && std::isfinite(values.eta) &&
           std::isfinite(values.phi) && std::isfinite(values.pt);
  }

  void fillPropagatedTrackHistograms(PropagatedTrackValues const& values, bool isCA)
  {
    fillHistogram(HIST("propagated/pv/all/pt"), values.pt);
    fillHistogram(HIST("propagated/pv/all/eta"), values.eta);
    fillHistogram(HIST("propagated/pv/all/phi"), values.phi);
    fillHistogram(HIST("propagated/pv/all/x"), values.x);
    fillHistogram(HIST("propagated/pv/all/y"), values.y);
    fillHistogram(HIST("propagated/pv/all/xy"), values.x, values.y);
    fillHistogram(HIST("propagated/pv/all/dcaX"), values.dcaX);
    fillHistogram(HIST("propagated/pv/all/dcaY"), values.dcaY);
    fillHistogram(HIST("propagated/pv/all/dcaXY"), values.dcaXY);

    if (isCA) {
      fillHistogram(HIST("propagated/pv/ca/pt"), values.pt);
      fillHistogram(HIST("propagated/pv/ca/eta"), values.eta);
      fillHistogram(HIST("propagated/pv/ca/phi"), values.phi);
      fillHistogram(HIST("propagated/pv/ca/x"), values.x);
      fillHistogram(HIST("propagated/pv/ca/y"), values.y);
      fillHistogram(HIST("propagated/pv/ca/xy"), values.x, values.y);
      fillHistogram(HIST("propagated/pv/ca/dcaX"), values.dcaX);
      fillHistogram(HIST("propagated/pv/ca/dcaY"), values.dcaY);
      fillHistogram(HIST("propagated/pv/ca/dcaXY"), values.dcaXY);
      return;
    }

    fillHistogram(HIST("propagated/pv/ltf/pt"), values.pt);
    fillHistogram(HIST("propagated/pv/ltf/eta"), values.eta);
    fillHistogram(HIST("propagated/pv/ltf/phi"), values.phi);
    fillHistogram(HIST("propagated/pv/ltf/x"), values.x);
    fillHistogram(HIST("propagated/pv/ltf/y"), values.y);
    fillHistogram(HIST("propagated/pv/ltf/xy"), values.x, values.y);
    fillHistogram(HIST("propagated/pv/ltf/dcaX"), values.dcaX);
    fillHistogram(HIST("propagated/pv/ltf/dcaY"), values.dcaY);
    fillHistogram(HIST("propagated/pv/ltf/dcaXY"), values.dcaXY);
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track) const
  {
    const float normalizedChi2 = chi2Ndf(track);
    return std::isfinite(track.pt()) && std::isfinite(track.eta()) && normalizedChi2 >= 0.f &&
           track.pt() >= selection.ptMin.value && track.pt() <= selection.ptMax.value &&
           track.eta() > selection.etaMin.value && track.eta() < selection.etaMax.value &&
           track.nClusters() >= selection.minNClusters.value && normalizedChi2 <= selection.maxChi2Ndf.value;
  }

  template <typename TParticle>
  bool isSelectedTruthParticle(TParticle const& particle) const
  {
    if (!std::isfinite(particle.pt()) || !std::isfinite(particle.eta()) ||
        particle.pt() < selection.ptMin.value || particle.pt() > selection.ptMax.value ||
        particle.eta() <= selection.etaMin.value || particle.eta() >= selection.etaMax.value) {
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
    fillHistogram(HIST("reco/all/pt"), track.pt());
    fillHistogram(HIST("reco/all/eta"), track.eta());
    fillHistogram(HIST("reco/all/phi"), phi);
    fillHistogram(HIST("reco/all/x"), track.x());
    fillHistogram(HIST("reco/all/y"), track.y());
    fillHistogram(HIST("reco/all/z"), track.z());
    fillHistogram(HIST("reco/all/signed1Pt"), track.signed1Pt());
    fillHistogram(HIST("reco/all/sign"), track.sign());
    fillHistogram(HIST("reco/all/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      fillHistogram(HIST("reco/all/chi2Ndf"), normalizedChi2);
      fillHistogram(HIST("reco/all/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    fillHistogram(HIST("reco/all/xy"), track.x(), track.y());
    fillHistogram(HIST("reco/all/etaPhi"), track.eta(), phi);
    fillHistogram(HIST("reco/all/nClustersVsEta"), track.eta(), track.nClusters());
    fillHistogram(HIST("reco/all/trackTime"), track.trackTime());
    fillHistogram(HIST("reco/all/trackTimeRes"), track.trackTimeRes());
  }

  template <typename TTrack>
  void fillAlgorithmTrackHistograms(TTrack const& track, float normalizedChi2, float phi)
  {
    if (track.isCA()) {
      fillHistogram(HIST("reco/ca/pt"), track.pt());
      fillHistogram(HIST("reco/ca/eta"), track.eta());
      fillHistogram(HIST("reco/ca/phi"), phi);
      fillHistogram(HIST("reco/ca/x"), track.x());
      fillHistogram(HIST("reco/ca/y"), track.y());
      fillHistogram(HIST("reco/ca/z"), track.z());
      fillHistogram(HIST("reco/ca/signed1Pt"), track.signed1Pt());
      fillHistogram(HIST("reco/ca/sign"), track.sign());
      fillHistogram(HIST("reco/ca/nClusters"), track.nClusters());
      if (normalizedChi2 >= 0.f) {
        fillHistogram(HIST("reco/ca/chi2Ndf"), normalizedChi2);
        fillHistogram(HIST("reco/ca/chi2NdfVsPt"), track.pt(), normalizedChi2);
      }
      fillHistogram(HIST("reco/ca/xy"), track.x(), track.y());
      fillHistogram(HIST("reco/ca/etaPhi"), track.eta(), phi);
      fillHistogram(HIST("reco/ca/nClustersVsEta"), track.eta(), track.nClusters());
      fillHistogram(HIST("reco/ca/trackTime"), track.trackTime());
      fillHistogram(HIST("reco/ca/trackTimeRes"), track.trackTimeRes());
      return;
    }

    fillHistogram(HIST("reco/ltf/pt"), track.pt());
    fillHistogram(HIST("reco/ltf/eta"), track.eta());
    fillHistogram(HIST("reco/ltf/phi"), phi);
    fillHistogram(HIST("reco/ltf/x"), track.x());
    fillHistogram(HIST("reco/ltf/y"), track.y());
    fillHistogram(HIST("reco/ltf/z"), track.z());
    fillHistogram(HIST("reco/ltf/signed1Pt"), track.signed1Pt());
    fillHistogram(HIST("reco/ltf/sign"), track.sign());
    fillHistogram(HIST("reco/ltf/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      fillHistogram(HIST("reco/ltf/chi2Ndf"), normalizedChi2);
      fillHistogram(HIST("reco/ltf/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    fillHistogram(HIST("reco/ltf/xy"), track.x(), track.y());
    fillHistogram(HIST("reco/ltf/etaPhi"), track.eta(), phi);
    fillHistogram(HIST("reco/ltf/nClustersVsEta"), track.eta(), track.nClusters());
    fillHistogram(HIST("reco/ltf/trackTime"), track.trackTime());
    fillHistogram(HIST("reco/ltf/trackTimeRes"), track.trackTimeRes());
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
      fillHistogram(HIST("reco/all/layerOccupancy"), layer);
      fillHistogram(HIST("reco/all/clusterSizeVsLayer"), layer, clusterSize);
      if (track.isCA()) {
        fillHistogram(HIST("reco/ca/layerOccupancy"), layer);
        fillHistogram(HIST("reco/ca/clusterSizeVsLayer"), layer, clusterSize);
      } else {
        fillHistogram(HIST("reco/ltf/layerOccupancy"), layer);
        fillHistogram(HIST("reco/ltf/clusterSizeVsLayer"), layer, clusterSize);
      }
    }
    fillHistogram(HIST("reco/all/layerPattern"), layerPattern);
    if (track.isCA()) {
      fillHistogram(HIST("reco/ca/layerPattern"), layerPattern);
    } else {
      fillHistogram(HIST("reco/ltf/layerPattern"), layerPattern);
    }
  }

  template <typename TTrack>
  void fillSelectedTrackHistograms(TTrack const& track, float normalizedChi2, float phi)
  {
    fillHistogram(HIST("reco/selected/all/pt"), track.pt());
    fillHistogram(HIST("reco/selected/all/eta"), track.eta());
    fillHistogram(HIST("reco/selected/all/phi"), phi);
    fillHistogram(HIST("reco/selected/all/x"), track.x());
    fillHistogram(HIST("reco/selected/all/y"), track.y());
    fillHistogram(HIST("reco/selected/all/z"), track.z());
    fillHistogram(HIST("reco/selected/all/signed1Pt"), track.signed1Pt());
    fillHistogram(HIST("reco/selected/all/sign"), track.sign());
    fillHistogram(HIST("reco/selected/all/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      fillHistogram(HIST("reco/selected/all/chi2Ndf"), normalizedChi2);
      fillHistogram(HIST("reco/selected/all/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    fillHistogram(HIST("reco/selected/all/xy"), track.x(), track.y());
    fillHistogram(HIST("reco/selected/all/etaPhi"), track.eta(), phi);
    fillHistogram(HIST("reco/selected/all/nClustersVsEta"), track.eta(), track.nClusters());
    fillHistogram(HIST("reco/selected/all/trackTime"), track.trackTime());
    fillHistogram(HIST("reco/selected/all/trackTimeRes"), track.trackTimeRes());

    if (track.isCA()) {
      fillHistogram(HIST("reco/selected/ca/pt"), track.pt());
      fillHistogram(HIST("reco/selected/ca/eta"), track.eta());
      fillHistogram(HIST("reco/selected/ca/phi"), phi);
      fillHistogram(HIST("reco/selected/ca/x"), track.x());
      fillHistogram(HIST("reco/selected/ca/y"), track.y());
      fillHistogram(HIST("reco/selected/ca/z"), track.z());
      fillHistogram(HIST("reco/selected/ca/signed1Pt"), track.signed1Pt());
      fillHistogram(HIST("reco/selected/ca/sign"), track.sign());
      fillHistogram(HIST("reco/selected/ca/nClusters"), track.nClusters());
      if (normalizedChi2 >= 0.f) {
        fillHistogram(HIST("reco/selected/ca/chi2Ndf"), normalizedChi2);
        fillHistogram(HIST("reco/selected/ca/chi2NdfVsPt"), track.pt(), normalizedChi2);
      }
      fillHistogram(HIST("reco/selected/ca/xy"), track.x(), track.y());
      fillHistogram(HIST("reco/selected/ca/etaPhi"), track.eta(), phi);
      fillHistogram(HIST("reco/selected/ca/nClustersVsEta"), track.eta(), track.nClusters());
      fillHistogram(HIST("reco/selected/ca/trackTime"), track.trackTime());
      fillHistogram(HIST("reco/selected/ca/trackTimeRes"), track.trackTimeRes());
      return;
    }

    fillHistogram(HIST("reco/selected/ltf/pt"), track.pt());
    fillHistogram(HIST("reco/selected/ltf/eta"), track.eta());
    fillHistogram(HIST("reco/selected/ltf/phi"), phi);
    fillHistogram(HIST("reco/selected/ltf/x"), track.x());
    fillHistogram(HIST("reco/selected/ltf/y"), track.y());
    fillHistogram(HIST("reco/selected/ltf/z"), track.z());
    fillHistogram(HIST("reco/selected/ltf/signed1Pt"), track.signed1Pt());
    fillHistogram(HIST("reco/selected/ltf/sign"), track.sign());
    fillHistogram(HIST("reco/selected/ltf/nClusters"), track.nClusters());
    if (normalizedChi2 >= 0.f) {
      fillHistogram(HIST("reco/selected/ltf/chi2Ndf"), normalizedChi2);
      fillHistogram(HIST("reco/selected/ltf/chi2NdfVsPt"), track.pt(), normalizedChi2);
    }
    fillHistogram(HIST("reco/selected/ltf/xy"), track.x(), track.y());
    fillHistogram(HIST("reco/selected/ltf/etaPhi"), track.eta(), phi);
    fillHistogram(HIST("reco/selected/ltf/nClustersVsEta"), track.eta(), track.nClusters());
    fillHistogram(HIST("reco/selected/ltf/trackTime"), track.trackTime());
    fillHistogram(HIST("reco/selected/ltf/trackTimeRes"), track.trackTimeRes());
  }

  template <typename TTrack>
  void fillSelectedClusterLayerHistograms(TTrack const& track)
  {
    const uint64_t clusterSizesAndFlags = track.mftClusterSizesAndTrackFlags();
    uint16_t layerPattern = 0;
    for (int layer = 0; layer < 10; ++layer) {
      const int clusterSize = static_cast<int>((clusterSizesAndFlags >> (6 * layer)) & 0x3fULL);
      if (clusterSize == 0) {
        continue;
      }
      layerPattern |= static_cast<uint16_t>(1U << layer);
      fillHistogram(HIST("reco/selected/all/layerOccupancy"), layer);
      fillHistogram(HIST("reco/selected/all/clusterSizeVsLayer"), layer, clusterSize);
      if (track.isCA()) {
        fillHistogram(HIST("reco/selected/ca/layerOccupancy"), layer);
        fillHistogram(HIST("reco/selected/ca/clusterSizeVsLayer"), layer, clusterSize);
      } else {
        fillHistogram(HIST("reco/selected/ltf/layerOccupancy"), layer);
        fillHistogram(HIST("reco/selected/ltf/clusterSizeVsLayer"), layer, clusterSize);
      }
    }
    fillHistogram(HIST("reco/selected/all/layerPattern"), layerPattern);
    if (track.isCA()) {
      fillHistogram(HIST("reco/selected/ca/layerPattern"), layerPattern);
    } else {
      fillHistogram(HIST("reco/selected/ltf/layerPattern"), layerPattern);
    }
  }

  AssociationCategory associationCategory(int nAssociations) const
  {
    if (nAssociations == 0) {
      return AssociationCategory::Unassigned;
    }
    return nAssociations == 1 ? AssociationCategory::Unique : AssociationCategory::Ambiguous;
  }

  HistogramRegistry& registryForAssociationCategory(AssociationCategory category)
  {
    if (category == AssociationCategory::Unique) {
      return registryAssociatedUnique;
    }
    if (category == AssociationCategory::Ambiguous) {
      return registryAssociatedAmbiguous;
    }
    return registryAssociatedUnassigned;
  }

  AssociationCounts countAssociations(aod::MFTTrackAssoc const& associations) const
  {
    AssociationCounts counts;
    counts.reserve(associations.size());
    for (const auto& association : associations) {
      ++counts[association.mfttrackId()];
    }
    return counts;
  }

  int getAssociationCount(AssociationCounts const& associationCounts, int64_t trackId) const
  {
    const auto found = associationCounts.find(trackId);
    return found == associationCounts.end() ? 0 : found->second;
  }

  template <typename TCollision>
  void updateBestCollision(TCollision const& collision,
                           bool propagationSucceeded,
                           PropagatedTrackValues const& propagatedValues,
                           BestCollisionResult& result)
  {
    const bool validCandidate = propagationSucceeded && std::isfinite(propagatedValues.dcaXY);
    fillHistogram(HIST("association/candidatePropagationStatus"), validCandidate ? 1. : 0.);
    if (!validCandidate) {
      return;
    }

    const int64_t collisionId = collision.globalIndex();
    const float dcaXY = propagatedValues.dcaXY;
    // Follow the MFT ambiguous-track reassociation convention: choose the
    // candidate with the smallest transverse distance after propagation.
    if (!result.found || dcaXY < result.dcaXY ||
        (dcaXY == result.dcaXY && collisionId < result.collisionId)) {
      result.found = true;
      result.collisionId = collisionId;
      result.dcaXY = dcaXY;
      result.propagatedValues = propagatedValues;
    }
  }

  template <typename TCollision>
  double collisionTimeDifference(TCollision const& otherCollision, TCollision const& bestCollision) const
  {
    // Positive values mean that the non-best collision occurred after the best collision.
    const int64_t bcOffset = static_cast<int64_t>(otherCollision.template bc_as<ExtBCs>().globalBC()) -
                             static_cast<int64_t>(bestCollision.template bc_as<ExtBCs>().globalBC());
    return static_cast<double>(bcOffset) * o2::constants::lhc::LHCBunchSpacingNS +
           static_cast<double>(otherCollision.collisionTime()) - static_cast<double>(bestCollision.collisionTime());
  }

  template <typename TCollision>
  void fillNonBestCollisionTimeDifference(TCollision const& otherCollision,
                                          TCollision const& bestCollision,
                                          bool isSelected)
  {
    const double deltaTime = collisionTimeDifference(otherCollision, bestCollision);
    fillHistogram(HIST("association/nonBestCollisionDeltaTime/all"), deltaTime);
    if (isSelected) {
      fillHistogram(HIST("association/nonBestCollisionDeltaTime/selected"), deltaTime);
    }
  }

  template <typename TTrack>
  bool fillTrackQA(TTrack const& track, bool isAssociated)
  {
    fillHistogram(HIST("reco/selectionCounter"), 0.);
    fillHistogram(HIST("reco/associationStatus"), isAssociated ? 1. : 0.);

    const float normalizedChi2 = chi2Ndf(track);
    const float phi = wrapPhi(track.phi());
    fillAllTrackHistograms(track, normalizedChi2, phi);
    fillAlgorithmTrackHistograms(track, normalizedChi2, phi);
    fillClusterLayerHistograms(track);

    if (!isSelectedTrack(track)) {
      return false;
    }
    fillHistogram(HIST("reco/selectionCounter"), 1.);
    fillHistogram(HIST("reco/selectedAlgorithm"), track.isCA() ? 1. : 0.);
    fillSelectedTrackHistograms(track, normalizedChi2, phi);
    fillSelectedClusterLayerHistograms(track);
    fillHistogram(HIST("caContribution/fractionVsPt"), track.pt(), track.isCA());
    fillHistogram(HIST("caContribution/fractionVsEta"), track.eta(), track.isCA());
    fillHistogram(HIST("caContribution/fractionVsNClusters"), track.nClusters(), track.isCA());
    fillHistogram(HIST("caContribution/fractionVsChi2Ndf"), normalizedChi2, track.isCA());
    return true;
  }

  template <typename TTrack, typename TCollision>
  bool propagateToCollisionPV(TTrack const& track, TCollision const& collision, PropagatedTrackValues& propagatedValues) const
  {
    bool propagationSucceeded = propagateToPV(track, collision.posZ(), propagatedValues);
    if (propagationSucceeded && std::isfinite(collision.posX()) && std::isfinite(collision.posY())) {
      propagatedValues.dcaX = propagatedValues.x - collision.posX();
      propagatedValues.dcaY = propagatedValues.y - collision.posY();
      propagatedValues.dcaXY = std::hypot(propagatedValues.dcaX, propagatedValues.dcaY);
      propagationSucceeded = std::isfinite(propagatedValues.dcaX) &&
                             std::isfinite(propagatedValues.dcaY) && std::isfinite(propagatedValues.dcaXY);
    } else {
      propagationSucceeded = false;
    }
    return propagationSucceeded;
  }

  template <typename TTrack, typename TCollision>
  bool fillPVPropagationQA(TTrack const& track, TCollision const& collision, PropagatedTrackValues& propagatedValues)
  {
    const bool propagationSucceeded = propagateToCollisionPV(track, collision, propagatedValues);
    fillHistogram(HIST("propagated/pv/status"), propagationSucceeded ? 2. : 1.);
    if (propagationSucceeded) {
      fillPropagatedTrackHistograms(propagatedValues, track.isCA());
    }
    return propagationSucceeded;
  }

  template <typename TTrack, typename TCollision>
  void fillPVPropagationQA(TTrack const& track, TCollision const& collision)
  {
    PropagatedTrackValues propagatedValues;
    fillPVPropagationQA(track, collision, propagatedValues);
  }

  void fillNoCollisionPropagationStatus()
  {
    fillHistogram(HIST("propagated/pv/status"), 0.);
  }

  template <typename TTrack>
  void addSelectedTrackToCollision(std::unordered_map<int64_t, TrackCounts>& countsPerCollision, int64_t collisionId, TTrack const& track)
  {
    auto& counts = countsPerCollision[collisionId];
    ++counts.all;
    if (track.isCA()) {
      ++counts.ca;
    } else {
      ++counts.ltf;
    }
  }

  template <typename TCollisions>
  void fillCollisionQA(TCollisions const& collisions, std::unordered_map<int64_t, TrackCounts> const& countsPerCollision)
  {
    for (const auto& collision : collisions) {
      const auto found = countsPerCollision.find(collision.globalIndex());
      const TrackCounts counts = found == countsPerCollision.end() ? TrackCounts{} : found->second;
      fillHistogram(HIST("collision/nTracksAll"), counts.all);
      fillHistogram(HIST("collision/nTracksLTF"), counts.ltf);
      fillHistogram(HIST("collision/nTracksCA"), counts.ca);
      fillHistogram(HIST("collision/nCAVsNLTF"), counts.ltf, counts.ca);
      if (counts.all > 0) {
        const float caFraction = static_cast<float>(counts.ca) / static_cast<float>(counts.all);
        fillHistogram(HIST("collision/caFraction"), caFraction);
        fillHistogram(HIST("collision/caFractionVsNTracks"), counts.all, caFraction);
      }
    }
  }

  template <typename TCollisions, typename TTracks>
  void fillRecoQA(TCollisions const& collisions, TTracks const& tracks)
  {
    std::unordered_map<int64_t, TrackCounts> countsPerCollision;
    for (const auto& track : tracks) {
      const bool hasCollision = track.has_collision();
      const bool isSelected = fillTrackQA(track, hasCollision);
      if (!hasCollision) {
        fillNoCollisionPropagationStatus();
        continue;
      }

      const auto collision = track.template collision_as<TCollisions>();
      fillPVPropagationQA(track, collision);
      if (isSelected) {
        addSelectedTrackToCollision(countsPerCollision, collision.globalIndex(), track);
      }
    }
    fillCollisionQA(collisions, countsPerCollision);
  }

  template <typename TCollisions, typename TTracks>
  BestCollisionResults fillAssociationCategorizedRecoQA(TCollisions const& collisions,
                                                        TTracks const& tracks,
                                                        aod::MFTTrackAssoc const& associations,
                                                        AssociationCounts const& associationCounts)
  {
    BestCollisionResults bestCollisions;
    bestCollisions.reserve(associationCounts.size());
    std::unordered_map<int64_t, TrackCounts> uniqueCountsPerCollision;
    std::unordered_map<int64_t, TrackCounts> ambiguousCountsPerCollision;
    for (const auto& association : associations) {
      const int nAssociations = getAssociationCount(associationCounts, association.mfttrackId());
      const auto category = associationCategory(nAssociations);
      activeRegistry = &registryForAssociationCategory(category);

      const auto track = association.template mfttrack_as<TTracks>();
      const auto collision = association.template collision_as<TCollisions>();
      PropagatedTrackValues propagatedValues;
      const bool propagationSucceeded = fillPVPropagationQA(track, collision, propagatedValues);
      if (category == AssociationCategory::Ambiguous) {
        activeRegistry = &registryAssociatedAmbiguousBest;
        updateBestCollision(collision, propagationSucceeded, propagatedValues, bestCollisions[association.mfttrackId()]);
      }

      if (!isSelectedTrack(track)) {
        continue;
      }
      auto& countsPerCollision = category == AssociationCategory::Unique ? uniqueCountsPerCollision : ambiguousCountsPerCollision;
      addSelectedTrackToCollision(countsPerCollision, collision.globalIndex(), track);
    }

    std::unordered_map<int64_t, TrackCounts> ambiguousBestCountsPerCollision;
    for (const auto& track : tracks) {
      const int nAssociations = getAssociationCount(associationCounts, track.globalIndex());
      const auto category = associationCategory(nAssociations);
      registryAssociationSummary.fill(HIST("nAssociationsPerTrack"), nAssociations);
      registryAssociationSummary.fill(HIST("trackCategory"), static_cast<uint8_t>(category));

      activeRegistry = &registryForAssociationCategory(category);
      fillTrackQA(track, nAssociations > 0);
      if (category == AssociationCategory::Unassigned) {
        fillNoCollisionPropagationStatus();
        continue;
      }
      if (category != AssociationCategory::Ambiguous) {
        continue;
      }

      activeRegistry = &registryAssociatedAmbiguousBest;
      const bool isSelected = fillTrackQA(track, true);
      const auto bestCollisionFound = bestCollisions.find(track.globalIndex());
      const bool hasBestCollision = bestCollisionFound != bestCollisions.end() && bestCollisionFound->second.found;
      fillHistogram(HIST("association/bestCollisionStatus"), hasBestCollision ? 1. : 0.);
      if (!hasBestCollision) {
        fillHistogram(HIST("propagated/pv/status"), 1.);
        continue;
      }

      const auto& bestCollision = bestCollisionFound->second;
      fillHistogram(HIST("association/bestDcaXY"), bestCollision.dcaXY);
      fillHistogram(HIST("propagated/pv/status"), 2.);
      fillPropagatedTrackHistograms(bestCollision.propagatedValues, track.isCA());
      if (isSelected) {
        addSelectedTrackToCollision(ambiguousBestCountsPerCollision, bestCollision.collisionId, track);
      }
    }

    for (const auto& association : associations) {
      const int nAssociations = getAssociationCount(associationCounts, association.mfttrackId());
      if (associationCategory(nAssociations) != AssociationCategory::Ambiguous) {
        continue;
      }
      const auto bestCollisionFound = bestCollisions.find(association.mfttrackId());
      if (bestCollisionFound == bestCollisions.end() || !bestCollisionFound->second.found ||
          association.collisionId() == bestCollisionFound->second.collisionId) {
        continue;
      }
      const auto track = association.template mfttrack_as<TTracks>();
      const auto otherCollision = association.template collision_as<TCollisions>();
      const auto bestCollision = collisions.iteratorAt(bestCollisionFound->second.collisionId);
      activeRegistry = &registryAssociatedAmbiguousBest;
      fillNonBestCollisionTimeDifference(otherCollision, bestCollision, isSelectedTrack(track));
    }

    activeRegistry = &registryAssociatedUnique;
    fillCollisionQA(collisions, uniqueCountsPerCollision);
    activeRegistry = &registryAssociatedAmbiguous;
    fillCollisionQA(collisions, ambiguousCountsPerCollision);
    activeRegistry = &registryAssociatedAmbiguousBest;
    fillCollisionQA(collisions, ambiguousBestCountsPerCollision);
    activeRegistry = &registryAssociatedUnassigned;
    const std::unordered_map<int64_t, TrackCounts> unassignedCountsPerCollision;
    fillCollisionQA(collisions, unassignedCountsPerCollision);
    activeRegistry = &registry;
    return bestCollisions;
  }

  ReconstructedCollisionCounts countReconstructedCollisionsByMcCollision(CollisionsMC const& collisions) const
  {
    ReconstructedCollisionCounts counts;
    counts.reserve(collisions.size());
    for (const auto& collision : collisions) {
      // Nreco is the multiplicity of reconstructed collisions carrying a valid
      // label to this MC collision. The collision mcMask is not used here.
      if (collision.has_mcCollision()) {
        ++counts[collision.mcCollisionId()];
      }
    }
    return counts;
  }

  int getReconstructedCollisionCount(ReconstructedCollisionCounts const& counts, int64_t mcCollisionId) const
  {
    const auto found = counts.find(mcCollisionId);
    return found == counts.end() ? 0 : found->second;
  }

  AssociationTruthStatus associationTruthStatus(MFTTracksMC::iterator const& track,
                                                 ReconstructedCollisionCounts const& reconstructedCollisionCounts,
                                                 bool associationCorrect) const
  {
    // Fake, noise and missing-label tracks cannot define a truth collision.
    if (labelStatus(track) != 0) {
      return AssociationTruthStatus::TrackTruthUnavailable;
    }

    const auto particle = track.template mcParticle_as<aod::McParticles>();
    const int nReconstructedCollisions = getReconstructedCollisionCount(reconstructedCollisionCounts, particle.mcCollisionId());
    if (nReconstructedCollisions == 0) {
      return AssociationTruthStatus::TruthCollisionNotReconstructed;
    }
    if (nReconstructedCollisions == 1) {
      return associationCorrect ? AssociationTruthStatus::TruthCollisionReconstructedAssociationCorrect : AssociationTruthStatus::TruthCollisionReconstructedAssociationWrong;
    }
    return associationCorrect ? AssociationTruthStatus::TruthCollisionSplitAssociationCorrect : AssociationTruthStatus::TruthCollisionSplitAssociationWrong;
  }

  void fillAssociationTruthStatus(MFTTracksMC::iterator const& track,
                                  AssociationTruthStatus status,
                                  bool hasDcaXY = false,
                                  float dcaXY = 0.f)
  {
    const auto statusValue = static_cast<uint8_t>(status);
    fillHistogram(HIST("mc/association/truthStatus/all"), statusValue);
    if (!isSelectedTrack(track)) {
      return;
    }
    fillHistogram(HIST("mc/association/truthStatus/selected"), statusValue);
    fillHistogram(HIST("mc/association/truthStatusVsPt"), track.pt(), statusValue);
    fillHistogram(HIST("mc/association/truthStatusVsEta"), track.eta(), statusValue);
    fillHistogram(HIST("mc/association/truthStatusVsPhi"), wrapPhi(track.phi()), statusValue);
    fillHistogram(HIST("mc/association/truthStatusVsChi2Ndf"), chi2Ndf(track), statusValue);
    if (hasDcaXY) {
      fillHistogram(HIST("mc/association/truthStatusVsDcaXY"), dcaXY, statusValue);
    }
  }

  void fillAssociationTruthQA(CollisionsMC const& collisions,
                              MFTTracksMC const& tracks,
                              aod::MFTTrackAssoc const& associations,
                              AssociationCounts const& associationCounts,
                              BestCollisionResults const& bestCollisions)
  {
    const auto reconstructedCollisionCounts = countReconstructedCollisionsByMcCollision(collisions);
    AssociationTruthSummaries associationTruthSummaries;
    associationTruthSummaries.reserve(associationCounts.size());
    for (const auto& association : associations) {
      const int nAssociations = getAssociationCount(associationCounts, association.mfttrackId());
      const auto category = associationCategory(nAssociations);
      const auto track = association.template mfttrack_as<MFTTracksMC>();
      const auto collision = association.template collision_as<CollisionsMC>();
      auto& summary = associationTruthSummaries[association.mfttrackId()];

      if (category == AssociationCategory::Unique && isSelectedTrack(track)) {
        PropagatedTrackValues propagatedValues;
        if (propagateToCollisionPV(track, collision, propagatedValues)) {
          summary.hasCandidateDcaXY = true;
          summary.candidateDcaXY = propagatedValues.dcaXY;
        }
      }

      if (labelStatus(track) == 0 && collision.has_mcCollision()) {
        const auto particle = track.template mcParticle_as<aod::McParticles>();
        if (collision.mcCollisionId() == particle.mcCollisionId()) {
          summary.hasTruthMatchedCandidate = true;
        }
      }
    }

    for (const auto& track : tracks) {
      const int nAssociations = getAssociationCount(associationCounts, track.globalIndex());
      const auto category = associationCategory(nAssociations);
      const auto summaryFound = associationTruthSummaries.find(track.globalIndex());
      const AssociationTruthSummary summary = summaryFound == associationTruthSummaries.end() ? AssociationTruthSummary{} : summaryFound->second;
      const auto categoryStatus = associationTruthStatus(track, reconstructedCollisionCounts, summary.hasTruthMatchedCandidate);
      activeRegistry = &registryForAssociationCategory(category);
      fillAssociationTruthStatus(track, categoryStatus,
                                 category == AssociationCategory::Unique && summary.hasCandidateDcaXY,
                                 summary.candidateDcaXY);

      if (category == AssociationCategory::Ambiguous) {
        const auto bestCollisionFound = bestCollisions.find(track.globalIndex());
        const bool hasBestCollision = bestCollisionFound != bestCollisions.end() && bestCollisionFound->second.found;
        bool bestAssociationCorrect = false;
        if (hasBestCollision && labelStatus(track) == 0) {
          const auto particle = track.template mcParticle_as<aod::McParticles>();
          const auto bestCollision = collisions.iteratorAt(bestCollisionFound->second.collisionId);
          bestAssociationCorrect = bestCollision.has_mcCollision() && bestCollision.mcCollisionId() == particle.mcCollisionId();
        }
        const auto bestStatus = associationTruthStatus(track, reconstructedCollisionCounts, bestAssociationCorrect);
        activeRegistry = &registryAssociatedAmbiguousBest;
        fillAssociationTruthStatus(track, bestStatus, hasBestCollision,
                                   hasBestCollision ? bestCollisionFound->second.dcaXY : 0.f);
      }
    }
    activeRegistry = &registry;
  }

  template <typename TParticle>
  void fillTruthKinematics(char category, TParticle const& particle)
  {
    switch (category) {
      case 'd':
        fillHistogram(HIST("mc/truth/denominator/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/denominator/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/denominator/ptEta"), particle.pt(), particle.eta());
        break;
      case 'i':
        fillHistogram(HIST("mc/truth/inclusive/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/inclusive/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/inclusive/ptEta"), particle.pt(), particle.eta());
        break;
      case 'l':
        fillHistogram(HIST("mc/truth/ltf/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/ltf/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/ltf/ptEta"), particle.pt(), particle.eta());
        break;
      case 'c':
        fillHistogram(HIST("mc/truth/ca/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/ca/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/ca/ptEta"), particle.pt(), particle.eta());
        break;
      case 'L':
        fillHistogram(HIST("mc/truth/ltfOnly/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/ltfOnly/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/ltfOnly/ptEta"), particle.pt(), particle.eta());
        break;
      case 'C':
        fillHistogram(HIST("mc/truth/caOnly/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/caOnly/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/caOnly/ptEta"), particle.pt(), particle.eta());
        break;
      case 'b':
        fillHistogram(HIST("mc/truth/both/pt"), particle.pt());
        fillHistogram(HIST("mc/truth/both/eta"), particle.eta());
        fillHistogram(HIST("mc/truth/both/ptEta"), particle.pt(), particle.eta());
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
        fillHistogram(HIST("mc/ca/ptRelativeResidualVsTruePt"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
      }
      fillHistogram(HIST("mc/ca/etaResidualVsTrueEta"), particle.eta(), track.eta() - particle.eta());
      fillHistogram(HIST("mc/ca/phiResidualVsTrueEta"), particle.eta(), phiResidual);
      fillHistogram(HIST("mc/ca/chargeCorrectVsTruePt"), particle.pt(), chargeCorrect);
      return;
    }

    if (particle.pt() > 0.f) {
      fillHistogram(HIST("mc/ltf/ptRelativeResidualVsTruePt"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
    }
    fillHistogram(HIST("mc/ltf/etaResidualVsTrueEta"), particle.eta(), track.eta() - particle.eta());
    fillHistogram(HIST("mc/ltf/phiResidualVsTrueEta"), particle.eta(), phiResidual);
    fillHistogram(HIST("mc/ltf/chargeCorrectVsTruePt"), particle.pt(), chargeCorrect);
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

  void fillMCQA(MFTTracksMC const& tracks,
                aod::McParticles const& particles,
                AssociationCounts const* associationCounts = nullptr,
                AssociationCategory requestedCategory = AssociationCategory::Unassigned)
  {
    std::unordered_map<int64_t, uint8_t> recoMaskByParticle;
    std::unordered_map<int64_t, std::array<int, 2>> recoCountsByParticle;

    for (const auto& track : tracks) {
      if (associationCounts != nullptr) {
        const int nAssociations = getAssociationCount(*associationCounts, track.globalIndex());
        if (associationCategory(nAssociations) != requestedCategory) {
          continue;
        }
      }
      if (!isSelectedTrack(track)) {
        continue;
      }

      const int status = labelStatus(track);
      const bool hasNonCorrectLabel = status != 0;
      fillHistogram(HIST("mc/labelStatusVsAlgorithm"), track.isCA() ? 1. : 0., status);
      if (track.isCA()) {
        fillHistogram(HIST("mc/ca/nonCorrectLabelFractionVsPt"), track.pt(), hasNonCorrectLabel);
        fillHistogram(HIST("mc/ca/nonCorrectLabelFractionVsEta"), track.eta(), hasNonCorrectLabel);
        fillHistogram(HIST("mc/ca/nonCorrectLabelFractionVsNClusters"), track.nClusters(), hasNonCorrectLabel);
      } else {
        fillHistogram(HIST("mc/ltf/nonCorrectLabelFractionVsPt"), track.pt(), hasNonCorrectLabel);
        fillHistogram(HIST("mc/ltf/nonCorrectLabelFractionVsEta"), track.eta(), hasNonCorrectLabel);
        fillHistogram(HIST("mc/ltf/nonCorrectLabelFractionVsNClusters"), track.nClusters(), hasNonCorrectLabel);
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
      fillHistogram(HIST("mc/truth/nRecoTracks"), recoCounts[0] + recoCounts[1]);
      fillHistogram(HIST("mc/truth/nRecoTracksLTF"), recoCounts[0]);
      fillHistogram(HIST("mc/truth/nRecoTracksCA"), recoCounts[1]);

      if (recoMask == 0) {
        fillHistogram(HIST("mc/truth/recoveryCategory"), 0.);
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
        fillHistogram(HIST("mc/truth/recoveryCategory"), 1.);
      } else if (recoMask == kCA) {
        fillTruthKinematics('C', particle);
        fillHistogram(HIST("mc/truth/recoveryCategory"), 2.);
      } else {
        fillTruthKinematics('b', particle);
        fillHistogram(HIST("mc/truth/recoveryCategory"), 3.);
      }
    }
  }

  void processData(aod::Collisions const& collisions, aod::MFTTracks const& tracks, ExtBCs const& bcs)
  {
    initMagneticField(bcs);
    activeRegistry = &registry;
    fillRecoQA(collisions, tracks);
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processData, "Process data MFT standalone-track QA", true);

  void processDataWithAssociations(aod::Collisions const& collisions,
                                   aod::MFTTracks const& tracks,
                                   aod::MFTTrackAssoc const& associations,
                                   ExtBCs const& bcs)
  {
    initMagneticField(bcs);
    activeRegistry = &registry;
    fillRecoQA(collisions, tracks);
    const auto associationCounts = countAssociations(associations);
    fillAssociationCategorizedRecoQA(collisions, tracks, associations, associationCounts);
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processDataWithAssociations, "Process data QA with MFTTrackAssoc categories", false);

  void processMC(aod::Collisions const& collisions, MFTTracksMC const& tracks, aod::McParticles const& particles, ExtBCs const& bcs)
  {
    initMagneticField(bcs);
    activeRegistry = &registry;
    fillRecoQA(collisions, tracks);
    fillMCQA(tracks, particles);
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processMC, "Process MC MFT standalone-track QA", false);

  void processMCWithAssociations(CollisionsMC const& collisions,
                                 MFTTracksMC const& tracks,
                                 aod::MFTTrackAssoc const& associations,
                                 aod::McParticles const& particles,
                                 ExtBCs const& bcs)
  {
    initMagneticField(bcs);
    activeRegistry = &registry;
    fillRecoQA(collisions, tracks);
    fillMCQA(tracks, particles);

    const auto associationCounts = countAssociations(associations);
    const auto bestCollisions = fillAssociationCategorizedRecoQA(collisions, tracks, associations, associationCounts);
    fillAssociationTruthQA(collisions, tracks, associations, associationCounts, bestCollisions);
    activeRegistry = &registryAssociatedUnique;
    fillMCQA(tracks, particles, &associationCounts, AssociationCategory::Unique);
    activeRegistry = &registryAssociatedAmbiguous;
    fillMCQA(tracks, particles, &associationCounts, AssociationCategory::Ambiguous);
    activeRegistry = &registryAssociatedAmbiguousBest;
    fillMCQA(tracks, particles, &associationCounts, AssociationCategory::Ambiguous);
    activeRegistry = &registryAssociatedUnassigned;
    fillMCQA(tracks, particles, &associationCounts, AssociationCategory::Unassigned);
    activeRegistry = &registry;
  }
  PROCESS_SWITCH(MFTStandaloneTrackQa, processMCWithAssociations, "Process MC QA with MFTTrackAssoc categories", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTStandaloneTrackQa>(cfgc)};
}
