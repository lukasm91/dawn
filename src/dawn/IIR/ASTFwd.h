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

#ifndef DAWN_IIR_ASTFWD_H
#define DAWN_IIR_ASTFWD_H

#include "dawn/AST/ASTFwd.h"

namespace dawn {
namespace iir {
//
// TODO refactor_AST: this is TEMPORARY, should be changed in the future to template specialization
//
using AST = ast::AST;

using Stmt = ast::Stmt;
using BlockStmt = ast::BlockStmt;
using ExprStmt = ast::ExprStmt;
using ReturnStmt = ast::ReturnStmt;
using VarDeclStmt = ast::VarDeclStmt;
using VerticalRegionDeclStmt = ast::VerticalRegionDeclStmt;
using StencilCallDeclStmt = ast::StencilCallDeclStmt;
using BoundaryConditionDeclStmt = ast::BoundaryConditionDeclStmt;
using IfStmt = ast::IfStmt;
using ReductionOverNeighborStmt = ast::ReductionOverNeighborStmt;

using Expr = ast::Expr;
using NOPExpr = ast::NOPExpr;
using UnaryOperator = ast::UnaryOperator;
using BinaryOperator = ast::BinaryOperator;
using AssignmentExpr = ast::AssignmentExpr;
using TernaryOperator = ast::TernaryOperator;
using FunCallExpr = ast::FunCallExpr;
using StencilFunCallExpr = ast::StencilFunCallExpr;
using StencilFunArgExpr = ast::StencilFunArgExpr;
using VarAccessExpr = ast::VarAccessExpr;
using FieldAccessExpr = ast::FieldAccessExpr;
using LiteralAccessExpr = ast::LiteralAccessExpr;

using ASTHelper = ast::ASTHelper;
using ASTVisitor = ast::ASTVisitor;
} // namespace iir
} // namespace dawn

#endif
