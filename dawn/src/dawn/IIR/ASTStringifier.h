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

#ifndef DAWN_IIR_ASTSTRINGIFER_H
#define DAWN_IIR_ASTSTRINGIFER_H

#include "dawn/AST/ASTStringifier.h"
#include "dawn/IIR/ASTFwd.h"

namespace dawn {
namespace iir {
//
// TODO refactor_AST: this is TEMPORARY, will be removed in the future
//
using ASTStringifier = ast::ASTStringifier;
inline std::ostream& operator<<(std::ostream& os, const AST& ast) {
  return ast::operator<<(os, ast);
}
inline std::ostream& operator<<(std::ostream& os, const std::shared_ptr<Stmt>& expr) {
  return ast::operator<<(os, expr);
}
inline std::ostream& operator<<(std::ostream& os, const std::shared_ptr<Expr>& stmt) {
  return ast::operator<<(os, stmt);
}
} // namespace iir
} // namespace dawn

#endif
