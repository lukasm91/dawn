//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "dawn/IIR/Extents.h"
#include <gtest/gtest.h>

using namespace dawn;
using namespace iir;
namespace {

static void compareExtents(iir::Extents extents, const std::array<int, 6>& ref) {
  const auto& hExtents =
      dawn::iir::extent_cast<dawn::iir::CartesianExtent const&>(extents.horizontalExtent());
  const auto& vExtents = extents.verticalExtent();

  EXPECT_EQ(hExtents.iMinus(), ref[0]);
  EXPECT_EQ(hExtents.iPlus(), ref[1]);
  EXPECT_EQ(hExtents.jMinus(), ref[2]);
  EXPECT_EQ(hExtents.jPlus(), ref[3]);
  EXPECT_EQ(vExtents.minus(), ref[4]);
  EXPECT_EQ(vExtents.plus(), ref[5]);
}

TEST(ExtentsTest, Construction) {
  Extents extents(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, -1, 1, 0});
  compareExtents(extents, {-1, -1, 1, 1, 0, 0});
}

TEST(ExtentsTest, PointWise) {
  Extents extents1(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, 0, 1, 0});
  EXPECT_FALSE(extents1.isHorizontalPointwise());
  EXPECT_TRUE(extents1.isVerticalPointwise());

  Extents extents2(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, 0, 0, 1});
  EXPECT_TRUE(extents2.isHorizontalPointwise());
  EXPECT_FALSE(extents2.isVerticalPointwise());
}

TEST(ExtentsTest, Merge1) {
  Extents extents(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, -1, 1, 0});
  Extents extentsToMerge(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, 3, 2, 1});
  extents.merge(extentsToMerge);

  compareExtents(extents, {-1, 3, 1, 2, 0, 1});
}

TEST(ExtentsTest, Merge2) {
  Extents extents(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, -1, 1, 0});
  Extents extentsToMerge(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, -2, 2, 0});
  extents.merge(extentsToMerge);

  compareExtents(extents, {-2, -1, 1, 2, 0, 0});
}

TEST(ExtentsTest, Merge3) {
  Extents extents(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, -1, 1, 0});
  extents.merge(ast::Offsets{ast::cartesian, -2, 0, 0});

  compareExtents(extents, {-2, -1, 0, 1, 0, 0});
}

TEST(ExtentsTest, Add) {
  Extents extents(dawn::ast::cartesian_{}, -2, 2, 0, 0, 0, 0);
  auto addedExtents =
      dawn::iir::Extents::add(extents, Extents(dawn::ast::cartesian_{}, -2, 2, 0, 0, 0, 0));

  compareExtents(addedExtents, {-4, 4, 0, 0, 0, 0});
}

TEST(ExtentsTest, addCenter) {
  Extents extents(dawn::ast::cartesian_{}, 1, 1, -2, -2, 3, 3);
  extents.addVerticalCenter();
  compareExtents(extents, {1, 1, -2, -2, 0, 3});
}

TEST(ExtentsTest, Stringify) {
  Extents extents(dawn::ast::cartesian_{}, ast::Offsets{ast::cartesian, 1, -1, 2});
  std::stringstream ss;
  ss << extents;
  EXPECT_STREQ(ss.str().c_str(), "[(1, 1), (-1, -1), (2, 2)]");
}

TEST(ExtentsTest, verticalLoopOrder) {
  Extents extents(dawn::ast::cartesian_{}, 0, 0, 0, 0, -1, 2);
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(
                  iir::LoopOrderKind::LK_Forward,
                  Extents::VerticalLoopOrderDir::VL_CounterLoopOrder, false)) == (Extent{1, 2}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(
                  iir::LoopOrderKind::LK_Forward,
                  Extents::VerticalLoopOrderDir::VL_CounterLoopOrder, true)) == (Extent{0, 2}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(iir::LoopOrderKind::LK_Forward,
                                                  Extents::VerticalLoopOrderDir::VL_InLoopOrder,
                                                  false)) == (Extent{-1, -1}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(iir::LoopOrderKind::LK_Forward,
                                                  Extents::VerticalLoopOrderDir::VL_InLoopOrder,
                                                  true)) == (Extent{-1, 0}));

  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(
                  iir::LoopOrderKind::LK_Backward,
                  Extents::VerticalLoopOrderDir::VL_CounterLoopOrder, false)) == (Extent{-1, -1}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(
                  iir::LoopOrderKind::LK_Backward,
                  Extents::VerticalLoopOrderDir::VL_CounterLoopOrder, true)) == (Extent{-1, 0}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(iir::LoopOrderKind::LK_Backward,
                                                  Extents::VerticalLoopOrderDir::VL_InLoopOrder,
                                                  false)) == (Extent{1, 2}));
  EXPECT_TRUE((extents.getVerticalLoopOrderExtent(iir::LoopOrderKind::LK_Backward,
                                                  Extents::VerticalLoopOrderDir::VL_InLoopOrder,
                                                  true)) == (Extent{0, 2}));
}
} // anonymous namespace
