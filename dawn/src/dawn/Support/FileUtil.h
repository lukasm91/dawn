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

#ifndef DAWN_SUPPORT_FILEUTIL_H
#define DAWN_SUPPORT_FILEUTIL_H

#include "dawn/Support/StringRef.h"

namespace dawn {

/// @brief Extract the filename from `path`
///
/// This will only work on UNIX like platforms.
///
/// @ingroup support
StringRef getFilename(StringRef path);

/// @brief Extract the extension from `filename`
/// @ingroup support
StringRef getExtension(StringRef filename);

/// @brief Extract the filename without extension from `path`
/// @ingroup support
StringRef getFilenameWithoutExtension(StringRef path);

} // namespace dawn

#endif
