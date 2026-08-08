// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// This software is distributed under the terms of the GNU General Public License v3.

#ifndef O2PHYSICS_COMMON_TASKS_MFTMCHMATCHINGCATEGORIES_H_
#define O2PHYSICS_COMMON_TASKS_MFTMCHMATCHINGCATEGORIES_H_

#include <cstdint>

namespace o2::analysis::mftmchmatching
{
enum class McMatchCategory : int {
  Unknown = -1,
  Correct = 0,
  FalseOther = 1,
  FalseNoMchLabel = 2,
  FalseNoMftLabel = 3,
  FalseDifferentMcCollision = 4,
  MuonMother = 5
};

enum class McMatchStudyCategory : int {
  Unknown = -1,
  CorrectPrimaryMuon = 0,
  CorrectPiKDecayMuon = 1,
  CorrectOtherMuon = 2,
  FalseWithTrueCounterpart = 3,
  NoTrueMftCounterpart = 4
};

template <typename T>
constexpr int serialized(T value)
{
  return static_cast<int>(value);
}

constexpr bool isCorrectStudyCategory(McMatchStudyCategory category)
{
  return category == McMatchStudyCategory::CorrectPrimaryMuon ||
         category == McMatchStudyCategory::CorrectPiKDecayMuon ||
         category == McMatchStudyCategory::CorrectOtherMuon;
}

constexpr bool isFalseStudyCategory(McMatchStudyCategory category)
{
  return category == McMatchStudyCategory::FalseWithTrueCounterpart ||
         category == McMatchStudyCategory::NoTrueMftCounterpart;
}

constexpr bool isCorrectPair(bool isSameMcParticle, bool isPiKMotherMatch)
{
  return isSameMcParticle || isPiKMotherMatch;
}

constexpr McMatchCategory classifyPair(bool correct,
                                       bool hasMchMcParticle,
                                       bool hasMftMcParticle,
                                       int64_t mchMcCollisionId,
                                       int64_t mftMcCollisionId)
{
  if (correct) {
    return McMatchCategory::Correct;
  }
  if (!hasMchMcParticle) {
    return McMatchCategory::FalseNoMchLabel;
  }
  if (!hasMftMcParticle) {
    return McMatchCategory::FalseNoMftLabel;
  }
  if (mchMcCollisionId >= 0 && mftMcCollisionId >= 0 && mchMcCollisionId != mftMcCollisionId) {
    return McMatchCategory::FalseDifferentMcCollision;
  }
  return McMatchCategory::FalseOther;
}

constexpr McMatchStudyCategory classifyStudyPair(bool correct,
                                                 bool isMchMuon,
                                                 bool hasPiKMother,
                                                 bool isPhysicalPrimary,
                                                 bool hasAllowedMftCounterpart)
{
  if (correct) {
    if (isMchMuon && hasPiKMother) {
      return McMatchStudyCategory::CorrectPiKDecayMuon;
    }
    if (isMchMuon && isPhysicalPrimary) {
      return McMatchStudyCategory::CorrectPrimaryMuon;
    }
    return McMatchStudyCategory::CorrectOtherMuon;
  }
  return hasAllowedMftCounterpart ? McMatchStudyCategory::FalseWithTrueCounterpart : McMatchStudyCategory::NoTrueMftCounterpart;
}

static_assert(serialized(McMatchCategory::Unknown) == -1);
static_assert(serialized(McMatchCategory::Correct) == 0);
static_assert(serialized(McMatchCategory::FalseOther) == 1);
static_assert(serialized(McMatchCategory::FalseNoMchLabel) == 2);
static_assert(serialized(McMatchCategory::FalseNoMftLabel) == 3);
static_assert(serialized(McMatchCategory::FalseDifferentMcCollision) == 4);
static_assert(serialized(McMatchCategory::MuonMother) == 5);
static_assert(serialized(McMatchStudyCategory::Unknown) == -1);
static_assert(serialized(McMatchStudyCategory::CorrectPrimaryMuon) == 0);
static_assert(serialized(McMatchStudyCategory::CorrectPiKDecayMuon) == 1);
static_assert(serialized(McMatchStudyCategory::CorrectOtherMuon) == 2);
static_assert(serialized(McMatchStudyCategory::FalseWithTrueCounterpart) == 3);
static_assert(serialized(McMatchStudyCategory::NoTrueMftCounterpart) == 4);

} // namespace o2::analysis::mftmchmatching

#endif // O2PHYSICS_COMMON_TASKS_MFTMCHMATCHINGCATEGORIES_H_
