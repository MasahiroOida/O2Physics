// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef O2PHYSICS_COMMON_TASKS_GLOBALTRACKMATCHINGSTUDYDATAMODEL_H_
#define O2PHYSICS_COMMON_TASKS_GLOBALTRACKMATCHINGSTUDYDATAMODEL_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod::globaltrackmatchingstudy
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(GlobalTrackId, globalTrackId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(MatchingPlaneZ, matchingPlaneZ, float);
DECLARE_SOA_COLUMN(ZMCH, zMCH, float);
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TglMCH, tglMCH, float);
DECLARE_SOA_COLUMN(Signed1PtMCH, signed1PtMCH, float);
DECLARE_SOA_COLUMN(PtMCH, ptMCH, float);
DECLARE_SOA_COLUMN(EtaMCH, etaMCH, float);
DECLARE_SOA_COLUMN(ZMFT, zMFT, float);
DECLARE_SOA_COLUMN(XMFT, xMFT, float);
DECLARE_SOA_COLUMN(YMFT, yMFT, float);
DECLARE_SOA_COLUMN(PhiMFT, phiMFT, float);
DECLARE_SOA_COLUMN(TglMFT, tglMFT, float);
DECLARE_SOA_COLUMN(Signed1PtMFT, signed1PtMFT, float);
DECLARE_SOA_COLUMN(PtMFT, ptMFT, float);
DECLARE_SOA_COLUMN(EtaMFT, etaMFT, float);
DECLARE_SOA_COLUMN(ZGlobal, zGlobal, float);
DECLARE_SOA_COLUMN(XGlobal, xGlobal, float);
DECLARE_SOA_COLUMN(YGlobal, yGlobal, float);
DECLARE_SOA_COLUMN(PhiGlobal, phiGlobal, float);
DECLARE_SOA_COLUMN(TglGlobal, tglGlobal, float);
DECLARE_SOA_COLUMN(Signed1PtGlobal, signed1PtGlobal, float);
DECLARE_SOA_COLUMN(PtGlobal, ptGlobal, float);
DECLARE_SOA_COLUMN(EtaGlobal, etaGlobal, float);
DECLARE_SOA_COLUMN(PtTruthMCH, ptTruthMCH, float);
DECLARE_SOA_COLUMN(EtaTruthMCH, etaTruthMCH, float);
DECLARE_SOA_COLUMN(PhiTruthMCH, phiTruthMCH, float);
DECLARE_SOA_COLUMN(PtTruthMFT, ptTruthMFT, float);
DECLARE_SOA_COLUMN(EtaTruthMFT, etaTruthMFT, float);
DECLARE_SOA_COLUMN(PhiTruthMFT, phiTruthMFT, float);
DECLARE_SOA_COLUMN(PtTruthGlobal, ptTruthGlobal, float);
DECLARE_SOA_COLUMN(EtaTruthGlobal, etaTruthGlobal, float);
DECLARE_SOA_COLUMN(PhiTruthGlobal, phiTruthGlobal, float);
DECLARE_SOA_COLUMN(DeltaX, deltaX, float);
DECLARE_SOA_COLUMN(DeltaY, deltaY, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaTgl, deltaTgl, float);
DECLARE_SOA_COLUMN(DeltaSigned1Pt, deltaSigned1Pt, float);
DECLARE_SOA_COLUMN(MatchingChi2, matchingChi2, float);
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, int);
DECLARE_SOA_COLUMN(McMatchStudyCategory, mcMatchStudyCategory, int);
DECLARE_SOA_COLUMN(McPdgCodeMCH, mcPdgCodeMCH, int);
DECLARE_SOA_COLUMN(McPdgCodeMFT, mcPdgCodeMFT, int);
DECLARE_SOA_COLUMN(McFirstMotherPdgCodeMCH, mcFirstMotherPdgCodeMCH, int);
} // namespace o2::aod::globaltrackmatchingstudy

namespace o2::aod
{
DECLARE_SOA_TABLE(GlobalTrackMatchingStudies, "AOD", "GTRKMATCHSTUDY",
                  globaltrackmatchingstudy::CollisionId,
                  globaltrackmatchingstudy::GlobalTrackId,
                  globaltrackmatchingstudy::MchTrackId,
                  globaltrackmatchingstudy::MftTrackId,
                  globaltrackmatchingstudy::MatchingPlaneZ,
                  globaltrackmatchingstudy::ZMCH,
                  globaltrackmatchingstudy::XMCH,
                  globaltrackmatchingstudy::YMCH,
                  globaltrackmatchingstudy::PhiMCH,
                  globaltrackmatchingstudy::TglMCH,
                  globaltrackmatchingstudy::Signed1PtMCH,
                  globaltrackmatchingstudy::PtMCH,
                  globaltrackmatchingstudy::EtaMCH,
                  globaltrackmatchingstudy::ZMFT,
                  globaltrackmatchingstudy::XMFT,
                  globaltrackmatchingstudy::YMFT,
                  globaltrackmatchingstudy::PhiMFT,
                  globaltrackmatchingstudy::TglMFT,
                  globaltrackmatchingstudy::Signed1PtMFT,
                  globaltrackmatchingstudy::PtMFT,
                  globaltrackmatchingstudy::EtaMFT,
                  globaltrackmatchingstudy::ZGlobal,
                  globaltrackmatchingstudy::XGlobal,
                  globaltrackmatchingstudy::YGlobal,
                  globaltrackmatchingstudy::PhiGlobal,
                  globaltrackmatchingstudy::TglGlobal,
                  globaltrackmatchingstudy::Signed1PtGlobal,
                  globaltrackmatchingstudy::PtGlobal,
                  globaltrackmatchingstudy::EtaGlobal,
                  globaltrackmatchingstudy::PtTruthMCH,
                  globaltrackmatchingstudy::EtaTruthMCH,
                  globaltrackmatchingstudy::PhiTruthMCH,
                  globaltrackmatchingstudy::PtTruthMFT,
                  globaltrackmatchingstudy::EtaTruthMFT,
                  globaltrackmatchingstudy::PhiTruthMFT,
                  globaltrackmatchingstudy::PtTruthGlobal,
                  globaltrackmatchingstudy::EtaTruthGlobal,
                  globaltrackmatchingstudy::PhiTruthGlobal,
                  globaltrackmatchingstudy::DeltaX,
                  globaltrackmatchingstudy::DeltaY,
                  globaltrackmatchingstudy::DeltaPhi,
                  globaltrackmatchingstudy::DeltaTgl,
                  globaltrackmatchingstudy::DeltaSigned1Pt,
                  globaltrackmatchingstudy::MatchingChi2,
                  globaltrackmatchingstudy::IsCorrectMatch,
                  globaltrackmatchingstudy::McMatchStudyCategory,
                  globaltrackmatchingstudy::McPdgCodeMCH,
                  globaltrackmatchingstudy::McPdgCodeMFT,
                  globaltrackmatchingstudy::McFirstMotherPdgCodeMCH);
} // namespace o2::aod

#endif // O2PHYSICS_COMMON_TASKS_GLOBALTRACKMATCHINGSTUDYDATAMODEL_H_
