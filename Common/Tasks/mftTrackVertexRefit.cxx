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

/// \file mftTrackVertexRefit.cxx
/// \brief Recompute MFT standalone track parameters with a collision vertex constraint.

#include "Common/DataModel/CollisionAssociationTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TDatabasePDG.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace mftvtxrefit
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);     //! Collision used as vertex constraint
DECLARE_SOA_INDEX_COLUMN(MFTTrack, mfttrack);       //! Original MFT track
DECLARE_SOA_COLUMN(Success, success, bool);         //! Vertex update status
DECLARE_SOA_COLUMN(X, x, float);                    //! Recomputed x at vertex z
DECLARE_SOA_COLUMN(Y, y, float);                    //! Recomputed y at vertex z
DECLARE_SOA_COLUMN(Z, z, float);                    //! Recomputed z
DECLARE_SOA_COLUMN(Phi, phi, float);                //! Recomputed phi
DECLARE_SOA_COLUMN(Tgl, tgl, float);                //! Recomputed tan(lambda)
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);    //! Recomputed q/pt
DECLARE_SOA_COLUMN(Pt, pt, float);                  //! Recomputed pt
DECLARE_SOA_COLUMN(Eta, eta, float);                //! Recomputed eta
DECLARE_SOA_COLUMN(P, p, float);                    //! Recomputed total momentum
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);      //! Recomputed phi - original phi
DECLARE_SOA_COLUMN(DeltaTgl, deltaTgl, float);      //! Recomputed tgl - original tgl
DECLARE_SOA_COLUMN(DeltaSigned1Pt, deltaSigned1Pt, float); //! Recomputed q/pt - original q/pt
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, float);        //! Recomputed pt - original pt
DECLARE_SOA_COLUMN(XOverX0, xOverX0, float);        //! Material budget used for MCS
} // namespace mftvtxrefit

DECLARE_SOA_TABLE(MFTVtxRefits, "AOD", "MFTVTXREFIT",
                  o2::soa::Index<>,
                  mftvtxrefit::CollisionId,
                  mftvtxrefit::MFTTrackId,
                  mftvtxrefit::Success,
                  mftvtxrefit::X,
                  mftvtxrefit::Y,
                  mftvtxrefit::Z,
                  mftvtxrefit::Phi,
                  mftvtxrefit::Tgl,
                  mftvtxrefit::Signed1Pt,
                  mftvtxrefit::Pt,
                  mftvtxrefit::Eta,
                  mftvtxrefit::P,
                  mftvtxrefit::DeltaPhi,
                  mftvtxrefit::DeltaTgl,
                  mftvtxrefit::DeltaSigned1Pt,
                  mftvtxrefit::DeltaPt,
                  mftvtxrefit::XOverX0);
} // namespace o2::aod

struct MFTTrackVertexRefit {
  Produces<aod::MFTVtxRefits> mftVtxRefits;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "CCDB path of the aligned geometry"};
  Configurable<double> bzInput{"bz", -999., "Bz field in kG. Values <= -998 load the field from CCDB"};
  Configurable<bool> useMaterialCorrection{"useMaterialCorrection", false, "Use mean material budget from geometry for MCS"};
  Configurable<float> xOverX0Fallback{"xOverX0Fallback", 0.f, "Fallback x/X0 used when material correction is disabled"};
  Configurable<float> fallbackSigmaX{"fallbackSigmaX", 0.1f, "Fallback MFT track sigma x in cm, used when MFTTracksCov is unavailable"};
  Configurable<float> fallbackSigmaY{"fallbackSigmaY", 0.1f, "Fallback MFT track sigma y in cm, used when MFTTracksCov is unavailable"};
  Configurable<float> fallbackSigmaPhi{"fallbackSigmaPhi", 0.1f, "Fallback MFT track sigma phi in rad, used when MFTTracksCov is unavailable"};
  Configurable<float> fallbackSigmaTgl{"fallbackSigmaTgl", 0.5f, "Fallback MFT track sigma tan(lambda), used when MFTTracksCov is unavailable"};
  Configurable<float> fallbackSigmaSigned1Pt{"fallbackSigmaSigned1Pt", 10.f, "Fallback MFT track sigma q/pt in c/GeV, used when MFTTracksCov is unavailable"};
  Configurable<int> nBinsTruthPt{"nBinsTruthPt", 200, "Number of truth pT bins for resolution histograms"};
  Configurable<float> minTruthPt{"minTruthPt", 0.f, "Minimum truth pT for resolution histograms"};
  Configurable<float> maxTruthPt{"maxTruthPt", 20.f, "Maximum truth pT for resolution histograms"};
  Configurable<int> nBinsPtResolution{"nBinsPtResolution", 400, "Number of pT resolution bins"};
  Configurable<float> minPtResolution{"minPtResolution", -2.f, "Minimum pT resolution"};
  Configurable<float> maxPtResolution{"maxPtResolution", 2.f, "Maximum pT resolution"};
  Configurable<int> nBinsTglResolution{"nBinsTglResolution", 400, "Number of tan(lambda) resolution bins"};
  Configurable<float> minTglResolution{"minTglResolution", -2.f, "Minimum tan(lambda) resolution"};
  Configurable<float> maxTglResolution{"maxTglResolution", 2.f, "Maximum tan(lambda) resolution"};
  Configurable<int> nBinsPhiResolution{"nBinsPhiResolution", 400, "Number of phi resolution bins"};
  Configurable<float> minPhiResolution{"minPhiResolution", -2.f, "Minimum phi resolution"};
  Configurable<float> maxPhiResolution{"maxPhiResolution", 2.f, "Maximum phi resolution"};
  Configurable<int> nBinsSigned1PtResolution{"nBinsSigned1PtResolution", 400, "Number of q/pT resolution bins"};
  Configurable<float> minSigned1PtResolution{"minSigned1PtResolution", -2.f, "Minimum q/pT resolution"};
  Configurable<float> maxSigned1PtResolution{"maxSigned1PtResolution", 2.f, "Maximum q/pT resolution"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;
  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

  int runNumber = -1;
  float bz = 0.f;
  TDatabasePDG* pdg = nullptr;

  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdbUrl);
    pdg = TDatabasePDG::Instance();

    registry.add("hAssociationMultiplicity", "MFT association multiplicity;N_{coll};counts", kTH1F, {{10, -0.5, 9.5}});
    registry.add("hDeltaPt", "vertex-constrained #Delta p_{T};#Delta p_{T} (GeV/#it{c});counts", kTH1F, {{400, -2., 2.}});
    registry.add("hDeltaTgl", "vertex-constrained #Delta tan#lambda;#Delta tan#lambda;counts", kTH1F, {{400, -2., 2.}});

    AxisSpec axisTruthPt{nBinsTruthPt, minTruthPt, maxTruthPt, "p_{T}^{true} (GeV/#it{c})"};
    AxisSpec axisPtResolution{nBinsPtResolution, minPtResolution, maxPtResolution, "(p_{T}^{true} - p_{T}^{reco}) / p_{T}^{true}"};
    AxisSpec axisTglResolution{nBinsTglResolution, minTglResolution, maxTglResolution, "(tan#lambda^{true} - tan#lambda^{reco}) / tan#lambda^{true}"};
    AxisSpec axisPhiResolution{nBinsPhiResolution, minPhiResolution, maxPhiResolution, "(#varphi^{true} - #varphi^{reco}) / #varphi^{true}"};
    AxisSpec axisSigned1PtResolution{nBinsSigned1PtResolution, minSigned1PtResolution, maxSigned1PtResolution, "((q/p_{T})^{true} - (q/p_{T})^{reco}) / (q/p_{T})^{true}"};
    registry.add("hPtResolutionMFT", "MFT standalone p_{T} resolution;p_{T}^{true} (GeV/#it{c});(p_{T}^{true} - p_{T}^{reco}) / p_{T}^{true}", kTH2F, {axisTruthPt, axisPtResolution});
    registry.add("hPtResolutionMFTVertex", "vertex-constrained MFT p_{T} resolution;p_{T}^{true} (GeV/#it{c});(p_{T}^{true} - p_{T}^{reco}) / p_{T}^{true}", kTH2F, {axisTruthPt, axisPtResolution});
    registry.add("hTglResolutionMFT", "MFT standalone tan#lambda resolution;p_{T}^{true} (GeV/#it{c});(tan#lambda^{true} - tan#lambda^{reco}) / tan#lambda^{true}", kTH2F, {axisTruthPt, axisTglResolution});
    registry.add("hTglResolutionMFTVertex", "vertex-constrained MFT tan#lambda resolution;p_{T}^{true} (GeV/#it{c});(tan#lambda^{true} - tan#lambda^{reco}) / tan#lambda^{true}", kTH2F, {axisTruthPt, axisTglResolution});
    registry.add("hPhiResolutionMFT", "MFT standalone #varphi resolution;p_{T}^{true} (GeV/#it{c});(#varphi^{true} - #varphi^{reco}) / #varphi^{true}", kTH2F, {axisTruthPt, axisPhiResolution});
    registry.add("hPhiResolutionMFTVertex", "vertex-constrained MFT #varphi resolution;p_{T}^{true} (GeV/#it{c});(#varphi^{true} - #varphi^{reco}) / #varphi^{true}", kTH2F, {axisTruthPt, axisPhiResolution});
    registry.add("hSigned1PtResolutionMFT", "MFT standalone q/p_{T} resolution;p_{T}^{true} (GeV/#it{c});((q/p_{T})^{true} - (q/p_{T})^{reco}) / (q/p_{T})^{true}", kTH2F, {axisTruthPt, axisSigned1PtResolution});
    registry.add("hSigned1PtResolutionMFTVertex", "vertex-constrained MFT q/p_{T} resolution;p_{T}^{true} (GeV/#it{c});((q/p_{T})^{true} - (q/p_{T})^{reco}) / (q/p_{T})^{true}", kTH2F, {axisTruthPt, axisSigned1PtResolution});
  }

  void initFieldAndGeometry(int run)
  {
    if (run == runNumber) {
      return;
    }
    runNumber = run;

    if (bzInput > -998.) {
      bz = bzInput;
    } else {
      std::map<std::string, std::string> metadata;
      auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, runNumber);
      auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, soreor.first);
      o2::base::Propagator::initFieldFromGRP(grpmag);

      const double centerMFT[3] = {0., 0., -61.4};
      auto* field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      if (field == nullptr) {
        LOGF(fatal, "Magnetic field was not initialized");
      }
      bz = field->getBz(centerMFT);
    }

    if (useMaterialCorrection && !o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    LOGF(info, "MFT vertex refit: run=%d, Bz=%f kG", runNumber, bz);
  }

  static float getPt(float signed1Pt)
  {
    constexpr float almostZero = 1.e-12f;
    return std::abs(signed1Pt) < almostZero ? std::numeric_limits<float>::max() : std::abs(1.f / signed1Pt);
  }

  static float getEta(float tgl)
  {
    return -std::log(std::tan(o2::constants::math::PIQuarter - 0.5f * std::atan(tgl)));
  }

  static float getP(float pt, float tgl)
  {
    return pt * std::sqrt(1.f + tgl * tgl);
  }

  static float getPtResolution(float truthPt, float recoPt)
  {
    return (truthPt - recoPt) / truthPt;
  }

  static float getTglResolution(float truthTgl, float recoTgl)
  {
    return (truthTgl - recoTgl) / truthTgl;
  }

  static float getPhiResolution(float truthPhi, float recoPhi)
  {
    const float twoPi = 2.f * o2::constants::math::PI;
    return std::remainder(truthPhi - recoPhi, twoPi) / truthPhi;
  }

  static float getSigned1PtResolution(float truthSigned1Pt, float recoSigned1Pt)
  {
    return (truthSigned1Pt - recoSigned1Pt) / truthSigned1Pt;
  }

  template <typename TTrack, typename TTrackCov>
  static o2::track::TrackParCovFwd makeTrackParCovFwd(TTrack const& track, TTrackCov const& cov)
  {
    SMatrix5 pars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> covValues{
      cov.cXX(), cov.cXY(), cov.cYY(), cov.cPhiX(), cov.cPhiY(),
      cov.cPhiPhi(), cov.cTglX(), cov.cTglY(), cov.cTglPhi(), cov.cTglTgl(),
      cov.c1PtX(), cov.c1PtY(), cov.c1PtPhi(), cov.c1PtTgl(), cov.c1Pt21Pt2()};
    SMatrix55 covs(covValues.begin(), covValues.end());
    return o2::track::TrackParCovFwd(track.z(), pars, covs, track.chi2());
  }

  template <typename TTrack>
  o2::track::TrackParCovFwd makeTrackParCovFwd(TTrack const& track) const
  {
    SMatrix5 pars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    SMatrix55 covs;
    covs(0, 0) = fallbackSigmaX.value * fallbackSigmaX.value;
    covs(1, 1) = fallbackSigmaY.value * fallbackSigmaY.value;
    covs(2, 2) = fallbackSigmaPhi.value * fallbackSigmaPhi.value;
    covs(3, 3) = fallbackSigmaTgl.value * fallbackSigmaTgl.value;
    covs(4, 4) = fallbackSigmaSigned1Pt.value * fallbackSigmaSigned1Pt.value;
    return o2::track::TrackParCovFwd(track.z(), pars, covs, track.chi2());
  }

  template <typename TTrack, typename TCollision>
  float getXOverX0(TTrack const& track, TCollision const& collision) const
  {
    if (!useMaterialCorrection.value) {
      return xOverX0Fallback.value;
    }
    auto mat = o2::base::GeometryManager::meanMaterialBudget(track.x(), track.y(), track.z(), collision.posX(), collision.posY(), collision.posZ());
    return static_cast<float>(mat.meanX2X0);
  }

  template <typename TTrack, typename TCollision>
  void fillRefitTable(TTrack const& track, TCollision const& collision, o2::track::TrackParCovFwd& trackPar, float truthPt, float truthPhi, float truthTgl, float truthSigned1Pt)
  {
    const float originalPt = track.pt();
    trackPar.propagateToZhelix(collision.posZ(), bz);
    const bool success = trackPar.update({collision.posX(), collision.posY()},
                                         {collision.covXX(), collision.covYY()});

    const float refitSigned1Pt = trackPar.getInvQPt();
    const float refitPt = getPt(refitSigned1Pt);
    const float refitTgl = trackPar.getTgl();
    const float deltaPt = refitPt - originalPt;
    const float deltaTgl = refitTgl - track.tgl();

    registry.fill(HIST("hDeltaPt"), deltaPt);
    registry.fill(HIST("hDeltaTgl"), deltaTgl);
    if (truthPt > 0.f && std::isfinite(truthPt) && std::isfinite(originalPt)) {
      registry.fill(HIST("hPtResolutionMFT"), truthPt, getPtResolution(truthPt, originalPt));
      if (success && std::isfinite(refitPt)) {
        registry.fill(HIST("hPtResolutionMFTVertex"), truthPt, getPtResolution(truthPt, refitPt));
      }
    }
    if (truthPt > 0.f && std::isfinite(truthPt) && std::isfinite(truthTgl) && std::abs(truthTgl) > 1.e-12f && std::isfinite(track.tgl())) {
      registry.fill(HIST("hTglResolutionMFT"), truthPt, getTglResolution(truthTgl, track.tgl()));
      if (success && std::isfinite(refitTgl)) {
        registry.fill(HIST("hTglResolutionMFTVertex"), truthPt, getTglResolution(truthTgl, refitTgl));
      }
    }
    if (truthPt > 0.f && std::isfinite(truthPt) && std::isfinite(truthPhi) && std::abs(truthPhi) > 1.e-12f && std::isfinite(track.phi())) {
      registry.fill(HIST("hPhiResolutionMFT"), truthPt, getPhiResolution(truthPhi, track.phi()));
      if (success && std::isfinite(trackPar.getPhi())) {
        registry.fill(HIST("hPhiResolutionMFTVertex"), truthPt, getPhiResolution(truthPhi, trackPar.getPhi()));
      }
    }
    if (truthPt > 0.f && std::isfinite(truthPt) && std::isfinite(truthSigned1Pt) && std::abs(truthSigned1Pt) > 1.e-12f && std::isfinite(track.signed1Pt())) {
      registry.fill(HIST("hSigned1PtResolutionMFT"), truthPt, getSigned1PtResolution(truthSigned1Pt, track.signed1Pt()));
      if (success && std::isfinite(refitSigned1Pt)) {
        registry.fill(HIST("hSigned1PtResolutionMFTVertex"), truthPt, getSigned1PtResolution(truthSigned1Pt, refitSigned1Pt));
      }
    }

    mftVtxRefits(collision.globalIndex(),
                 track.globalIndex(),
                 success,
                 trackPar.getX(),
                 trackPar.getY(),
                 trackPar.getZ(),
                 trackPar.getPhi(),
                 refitTgl,
                 refitSigned1Pt,
                 refitPt,
                 getEta(refitTgl),
                 getP(refitPt, refitTgl),
                 trackPar.getPhi() - track.phi(),
                 deltaTgl,
                 refitSigned1Pt - track.signed1Pt(),
                 deltaPt,
                 0.f);
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

  void processNoCov(aod::Collisions const&,
                    aod::BCs const& bcs,
                    MFTTracksLabeled const&,
                    aod::MFTTrackAssoc const& associations,
                    aod::McParticles const&)
  {
    if (bcs.begin() == bcs.end()) {
      return;
    }
    initFieldAndGeometry(bcs.begin().runNumber());

    auto nAssociationsPerTrack = countAssociationsPerTrack(associations);

    mftVtxRefits.reserve(associations.size());
    for (auto const& association : associations) {
      auto nAssoc = nAssociationsPerTrack[association.mfttrackId()];
      registry.fill(HIST("hAssociationMultiplicity"), nAssoc);
      if (nAssoc != 1) {
        continue;
      }

      auto track = association.mfttrack_as<MFTTracksLabeled>();
      auto collision = association.collision_as<aod::Collisions>();
      float truthPt = -1.f;
      float truthPhi = std::numeric_limits<float>::quiet_NaN();
      float truthTgl = std::numeric_limits<float>::quiet_NaN();
      float truthSigned1Pt = std::numeric_limits<float>::quiet_NaN();
      if (track.has_mcParticle()) {
        auto mcParticle = track.template mcParticle_as<aod::McParticles>();
        truthPt = mcParticle.pt();
        if (truthPt > 0.f) {
          truthPhi = std::atan2(mcParticle.py(), mcParticle.px());
          truthTgl = mcParticle.pz() / truthPt;
          if (pdg != nullptr) {
            auto* pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
            if (pdgParticle != nullptr) {
              const float charge = pdgParticle->Charge() / 3.f;
              if (std::abs(charge) > 1.e-12f) {
                truthSigned1Pt = charge / truthPt;
              }
            }
          }
        }
      }
      auto trackPar = makeTrackParCovFwd(track);
      fillRefitTable(track, collision, trackPar, truthPt, truthPhi, truthTgl, truthSigned1Pt);
    }
  }
  PROCESS_SWITCH(MFTTrackVertexRefit, processNoCov, "Process MFT tracks without MFTTracksCov table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTTrackVertexRefit>(cfgc)};
}
