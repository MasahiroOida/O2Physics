# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2physics_add_dpl_workflow(trackqa
                    SOURCES trackqa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(propagatorqa
                    SOURCES propagatorQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2Physics::trackSelectionRequest
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(integrationtest
                    SOURCES integrationTest.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(integrationtestccdb
                    SOURCES integrationTestCCDB.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(validation
                    SOURCES validation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-qa
                    SOURCES multiplicityQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-derived-qa
                    SOURCES multiplicityDerivedQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-qa
                    SOURCES centralityQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-data-model
                  SOURCES checkDataModel.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-data-model-mc
                  SOURCES checkDataModelMC.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(orbit-range
                  SOURCES orbitRangeTask.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ft0-qa
                    SOURCES ft0Qa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-track-vertex-refit
                    SOURCES mftTrackVertexRefit.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB
                                          O2::DataFormatsParameters O2::ReconstructionDataFormats
                                          O2::DetectorsBase O2::Field
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-track-pair-mass
                    SOURCES mftTrackPairMass.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-standalone-track-qa
                    SOURCES mftStandaloneTrackQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB
                                          O2::DataFormatsParameters O2::ReconstructionDataFormats
                                          O2::DetectorsBase O2::Field
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-forward-muon-vz
                    SOURCES mcForwardMuonVz.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-mch-mid-direct-efficiency
                    SOURCES mftMchMidDirectEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(global-track-matching-study
                    SOURCES globalTrackMatchingStudy.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB
                                          O2::DataFormatsParameters O2::DetectorsBase
                                          O2::GlobalTracking O2::MCHTracking
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(global-track-dimuon-study
                    SOURCES globalTrackDimuonStudy.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon-mc-resolution
                    SOURCES dimuon_mc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase
                                          O2Physics::AnalysisCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-global-track
                    SOURCES check_globalTrack.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-cpv
                    SOURCES cpvQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvectors-correction
                    SOURCES qVectorsCorrection.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-study
                    SOURCES centralityStudy.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-study-pp
                    SOURCES centralityStudypp.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-test
                    SOURCES flowTest.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-qa
                    SOURCES qaMuon.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::Field O2::DetectorsBase O2::DetectorsCommonDataFormats O2::MathUtils O2::MCHTracking O2::DataFormatsMCH O2::GlobalTracking O2::MCHBase O2::MCHGeometryTransformer O2::CommonUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdc-table-reader
                    SOURCES zdcTableReader.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdc-extra-table-reader
                    SOURCES zdcExtraTableReader.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
