#include "dawn/CodeGen/StencilFunctionAsBCGenerator.h"
#include "dawn/IIR/ASTExpr.h"
#include "dawn/IIR/StencilInstantiation.h"

namespace dawn {
namespace codegen {

std::string
StencilFunctionAsBCGenerator::getName(const std::shared_ptr<iir::VarDeclStmt>& stmt) const {
  return metadata_.getFieldNameFromAccessID(iir::getAccessID(stmt));
}

std::string StencilFunctionAsBCGenerator::getName(const std::shared_ptr<iir::Expr>& expr) const {
  return metadata_.getFieldNameFromAccessID(iir::getAccessID(expr));
}

void StencilFunctionAsBCGenerator::visit(const std::shared_ptr<iir::FieldAccessExpr>& expr) {
  expr->getName();
  auto getArgumentIndex = [&](const std::string& name) {
    size_t pos =
        std::distance(function_->Args.begin(),
                      std::find_if(function_->Args.begin(), function_->Args.end(),
                                   [&](const std::shared_ptr<sir::StencilFunctionArg>& arg) {
                                     return arg->Name == name;
                                   }));

    DAWN_ASSERT_MSG(pos < function_->Args.size(), "");
    return pos;
  };
  ss_ << dawn::format("data_field_%i(%s)", getArgumentIndex(expr->getName()),
                      toString(expr->getOffset(), ", ", [&](std::string const& name, int offset) {
                        return name + "+" + std::to_string(offset);
                      }));
}

void StencilFunctionAsBCGenerator::visit(const std::shared_ptr<iir::VarAccessExpr>& expr) {
  if(metadata_.isAccessType(iir::FieldAccessType::FAT_GlobalVariable, iir::getAccessID(expr)))
    ss_ << "m_globals.";

  ss_ << getName(expr);

  if(expr->isArrayAccess()) {
    ss_ << "[";
    expr->getIndex()->accept(*this);
    ss_ << "]";
  }
}

void BCGenerator::generate(const std::shared_ptr<iir::BoundaryConditionDeclStmt>& stmt) {
  const auto& h_extents = dawn::iir::extent_cast<dawn::iir::CartesianExtent const&>(
      metadata_.getBoundaryConditionExtentsFromBCStmt(stmt).horizontalExtent());
  const auto& v_extents = metadata_.getBoundaryConditionExtentsFromBCStmt(stmt).verticalExtent();
  int haloIMinus = abs(h_extents.iMinus());
  int haloIPlus = abs(h_extents.iPlus());
  int haloJMinus = abs(h_extents.jMinus());
  int haloJPlus = abs(h_extents.jPlus());
  int haloKMinus = abs(v_extents.minus());
  int haloKPlus = abs(v_extents.plus());
  std::string fieldname = stmt->getFields()[0];

  // Set up the halos
  std::string halosetup = dawn::format(
      "gridtools::array< gridtools::halo_descriptor, 3 > halos;\n"
      "halos[0] =gridtools::halo_descriptor(%i, %i, "
      "%s.get_storage_info_ptr()->template begin<0>(),%s.get_storage_info_ptr()->template "
      "end<0>(), "
      "%s.get_storage_info_ptr()->template total_length<0>());\nhalos[1] = "
      "gridtools::halo_descriptor(%i, "
      "%i, "
      "%s.get_storage_info_ptr()->template begin<1>(),%s.get_storage_info_ptr()->template "
      "end<1>(), "
      "%s.get_storage_info_ptr()->template total_length<1>());\nhalos[2] = "
      "gridtools::halo_descriptor(%i, "
      "%i, "
      "%s.get_storage_info_ptr()->template begin<2>(),%s.get_storage_info_ptr()->template "
      "end<2>(), "
      "%s.get_storage_info_ptr()->template total_length<2>());\n",
      haloIMinus, haloIPlus, fieldname, fieldname, fieldname, haloJMinus, haloJPlus, fieldname,
      fieldname, fieldname, haloKMinus, haloKPlus, fieldname, fieldname, fieldname);
  std::string makeView = "";

  // Create the views for the fields
  for(int i = 0; i < stmt->getFields().size(); ++i) {
    auto fieldName = stmt->getFields()[i];
    makeView +=
        dawn::format("auto %s_view = GT_BACKEND_DECISION_viewmaker(%s);\n", fieldName, fieldName);
  }
  std::string bcapply = "GT_BACKEND_DECISION_bcapply<" + stmt->getFunctor() + " >(halos, " +
                        stmt->getFunctor() + "()).apply(";
  for(int i = 0; i < stmt->getFields().size(); ++i) {
    bcapply += stmt->getFields()[i] + "_view";
    if(i < stmt->getFields().size() - 1) {
      bcapply += ", ";
    }
  }
  bcapply += ");\n";

  ss_ << halosetup;
  ss_ << makeView;
  ss_ << bcapply;
}
} // namespace codegen
} // namespace dawn
