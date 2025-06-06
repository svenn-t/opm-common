/*
  Copyright 2019 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>

#include <opm/input/eclipse/Schedule/MSW/SegmentMatcher.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunction.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQSet.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDT.hpp>

#include <cassert>
#include <charconv>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <unordered_set>
#include <variant>
#include <vector>

#include <fmt/format.h>

namespace {

bool is_udq_blacklist(const std::string& keyword)
{
    static const auto udq_blacklistkw = std::unordered_set<std::string> {
        "SUMTHIN", "SUMMARY", "RUNSUM",
    };

    return udq_blacklistkw.find(keyword) != udq_blacklistkw.end();
}

bool is_udq(const std::string& keyword)
{
    // Does 'keyword' match one of the patterns
    //   AU*, BU*, CU*, FU*, GU*, RU*, SU*, or WU*?
    using sz_t = std::string::size_type;

    return (keyword.size() > sz_t{1})
        && (keyword[1] == 'U')
        && ! is_udq_blacklist(keyword)
        && (keyword.find_first_of("WGFCRBSA") == sz_t{0});
}

Opm::UDQVarType init_type(const Opm::UDQTokenType token_type)
{
    if ((token_type == Opm::UDQTokenType::number) ||
        Opm::UDQ::scalarFunc(token_type))
    {
        return Opm::UDQVarType::SCALAR;
    }

    return Opm::UDQVarType::NONE;
}

} // Anonymous namespace

namespace Opm {

UDQASTNode::UDQASTNode()
    : UDQASTNode(UDQTokenType::error)
{}

UDQASTNode::UDQASTNode(const UDQTokenType type_arg)
    : type(type_arg)
{
    if ((this->type == UDQTokenType::error) ||
        (this->type == UDQTokenType::binary_op_add) ||
        (this->type == UDQTokenType::binary_op_sub))
    {
        return;
    }

    throw std::invalid_argument {
        "Single argument AST node constructor available only "
        "for error and binary addition/subtraction tokens"
    };
}

UDQASTNode::UDQASTNode(double numeric_value)
    : var_type(init_type(UDQTokenType::number))
    , type    (UDQTokenType::number)
    , value   (numeric_value)
{}

UDQASTNode::UDQASTNode(const UDQTokenType                       type_arg,
                       const std::variant<std::string, double>& value_arg)
    : var_type(init_type(type_arg))
    , type    (type_arg)
    , value   (value_arg)
{}

UDQASTNode::UDQASTNode(const UDQTokenType                       type_arg,
                       const std::variant<std::string, double>& value_arg,
                       const UDQASTNode&                        left_arg)
    : UDQASTNode(type_arg, value_arg)
{
    if (UDQ::scalarFunc(type_arg)) {
        this->var_type = UDQVarType::SCALAR;
    }
    else {
        this->var_type = left_arg.var_type;
    }

    this->left = std::make_unique<UDQASTNode>(left_arg);
}

UDQASTNode::UDQASTNode(const UDQTokenType                       type_arg,
                       const std::variant<std::string, double>& value_arg,
                       const UDQASTNode&                        left_arg,
                       const UDQASTNode&                        right_arg)
    : var_type(init_type(type_arg))
    , type    (type_arg)
    , value   (value_arg)
{
    this->set_left(left_arg);
    this->set_right(right_arg);
}

UDQASTNode::UDQASTNode(const UDQTokenType                       type_arg,
                       const std::variant<std::string, double>& value_arg,
                       const std::vector<std::string>&          selector_arg)
    : var_type(init_type(type_arg))
    , type    (type_arg)
    , value   (value_arg)
    , selector(selector_arg)
{
    if (type_arg == UDQTokenType::ecl_expr) {
        this->var_type = UDQ::targetType(std::get<std::string>(this->value), this->selector);
    }

    if ((this->var_type == UDQVarType::CONNECTION_VAR) ||
        (this->var_type == UDQVarType::REGION_VAR) ||
        (this->var_type == UDQVarType::AQUIFER_VAR) ||
        (this->var_type == UDQVarType::BLOCK_VAR))
    {
        throw std::invalid_argument {
            fmt::format("UDQ variable type {} is not "
                        "currently supported for non-scalar uses",
                        UDQ::typeName(this->var_type))
        };
    }
}

UDQASTNode UDQASTNode::serializationTestObject()
{
    UDQASTNode result;
    result.var_type = UDQVarType::REGION_VAR;
    result.type = UDQTokenType::error;
    result.value = "test1";
    result.selector = {"test2"};
    result.sign = -1;

    UDQASTNode left = result;
    result.left = std::make_shared<UDQASTNode>(left);

    return result;
}

UDQSet
UDQASTNode::eval(const UDQVarType  target_type,
                 const UDQContext& context) const
{
    if (this->type == UDQTokenType::ecl_expr) {
        return this->sign * this->eval_expression(context);
    }

    if (UDQ::scalarFunc(this->type)) {
        return this->sign * this->eval_scalar_function(target_type, context);
    }

    if (UDQ::elementalUnaryFunc(this->type)) {
        return this->sign * this->eval_elemental_unary_function(target_type, context);
    }

    if (UDQ::binaryFunc(this->type)) {
        return this->sign * this->eval_binary_function(target_type, context);
    }

    if (this->type == UDQTokenType::number) {
        return this->sign * this->eval_number(target_type, context);
    }

    throw std::invalid_argument {
        "Should not be here ... this->type: " + std::to_string(static_cast<int>(this->type))
    };
}

bool UDQASTNode::valid() const
{
    return this->type != UDQTokenType::error;
}

std::set<UDQTokenType> UDQASTNode::func_tokens() const
{
    auto tokens = std::set<UDQTokenType>{};
    this->func_tokens(tokens);

    return tokens;
}

void UDQASTNode::update_type(const UDQASTNode& arg)
{
    this->var_type = (this->var_type == UDQVarType::NONE)
        ? arg.var_type
        : UDQ::coerce(this->var_type, arg.var_type);
}

void UDQASTNode::set_left(const UDQASTNode& arg)
{
    this->left = std::make_unique<UDQASTNode>(arg);
    this->update_type(arg);
}

void UDQASTNode::set_right(const UDQASTNode& arg)
{
    this->right = std::make_unique<UDQASTNode>(arg);
    this->update_type(arg);
}

void UDQASTNode::scale(double sign_factor)
{
    this->sign *= sign_factor;
}

UDQASTNode* UDQASTNode::get_left() const
{
    return this->left.get();
}

UDQASTNode* UDQASTNode::get_right() const
{
    return this->right.get();
}

bool UDQASTNode::operator==(const UDQASTNode& data) const
{
    if ((this->left && !data.left) ||
        (!this->left && data.left))
    {
        return false;
    }

    if (this->left && !(*this->left == *data.left)) {
        return false;
    }

    if ((this->right && !data.right) ||
        (!this->right && data.right))
    {
        return false;
    }

    if (this->right && !(*this->right == *data.right)) {
        return false;
    }

    return (type == data.type)
        && (var_type == data.var_type)
        && (value == data.value)
        && (selector == data.selector)
        ;
}

void UDQASTNode::required_summary(std::unordered_set<std::string>& summary_keys) const
{
    if ((this->type == UDQTokenType::ecl_expr) &&
        std::holds_alternative<std::string>(this->value))
    {
        if (const auto& keyword = std::get<std::string>(this->value);
            !is_udq(keyword))
        {
            summary_keys.insert(keyword);
        }
    }

    if (this->left) {
        this->left->required_summary(summary_keys);
    }

    if (this->right) {
        this->right->required_summary(summary_keys);
    }
}

void UDQASTNode::requiredObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    if ((this->type == UDQTokenType::ecl_expr) &&
        std::holds_alternative<std::string>(this->value))
    {
        this->populateRequiredObjects(objects);
    }

    if (this->left != nullptr) {
        this->left->requiredObjects(objects);
    }

    if (this->right != nullptr) {
        this->right->requiredObjects(objects);
    }
}

UDQSet
UDQASTNode::eval_expression(const UDQContext& context) const
{
    const auto& string_value = std::get<std::string>(this->value);
    const auto data_type = UDQ::targetType(string_value);

    if (data_type == UDQVarType::WELL_VAR) {
        return this->eval_well_expression(string_value, context);
    }

    if (data_type == UDQVarType::GROUP_VAR) {
        return this->eval_group_expression(string_value, context);
    }

    if (data_type == UDQVarType::SEGMENT_VAR) {
        return this->eval_segment_expression(string_value, context);
    }

    if (data_type == UDQVarType::REGION_VAR) {
        return this->eval_region_expression(string_value, context);
    }

    if (data_type == UDQVarType::FIELD_VAR) {
        return UDQSet::scalar(string_value, context.get(string_value));
    }

    if (data_type == UDQVarType::TABLE_LOOKUP) {
        const auto param_type = UDQ::targetType(this->selector[0]);
        return this->eval_table_lookup(param_type, string_value, context);
    }

    if (const auto scalar = context.get(string_value); scalar.has_value()) {
        return UDQSet::scalar(string_value, scalar.value());
    }

    throw std::logic_error {
        "Should not be here: var_type: '"
        + UDQ::typeName(data_type)
        + "' stringvalue: '"
        + string_value + '\''
    };
}

UDQSet
UDQASTNode::eval_well_expression(const std::string& string_value,
                                 const UDQContext&  context) const
{
    const auto& all_wells = context.wells();

    if (this->selector.empty()) {
        auto res = UDQSet::wells(string_value, all_wells);

        for (const auto& well : all_wells) {
            res.assign(well, context.get_well_var(well, string_value));
        }

        return res;
    }

    const auto& well_pattern = this->selector.front();

    if (well_pattern.find('*') == std::string::npos) {
        // The right hand side is a fully qualified well name without any
        // '*', in this case the right hand side evaluates to a *scalar* -
        // and that scalar value is distributed among all the wells in the
        // result set.
        return UDQSet::scalar(string_value, context.get_well_var(well_pattern, string_value));
    }
    else {
        // The right hand side is a set of wells.  The result set will be
        // updated for all wells in the right hand set, wells missing in the
        // right hand set will be undefined in the result set.
        auto res = UDQSet::wells(string_value, all_wells);
        for (const auto& wname : context.wells(well_pattern)) {
            res.assign(wname, context.get_well_var(wname, string_value));
        }

        return res;
    }
}

UDQSet
UDQASTNode::eval_group_expression(const std::string& string_value,
                                  const UDQContext&  context) const
{
    if (! this->selector.empty()) {
        const std::string& group_pattern = this->selector[0];
        if (group_pattern.find("*") == std::string::npos) {
            return UDQSet::scalar(string_value, context.get_group_var(group_pattern, string_value));
        }

        throw std::logic_error("Group names with wildcards is not yet supported");
    }

    const auto groups = context.nonFieldGroups();

    auto res = UDQSet::groups(string_value, groups);
    for (const auto& group : groups) {
        res.assign(group, context.get_group_var(group, string_value));
    }

    return res;
}

UDQSet
UDQASTNode::eval_segment_expression(const std::string& string_value,
                                    const UDQContext&  context) const
{
    const auto all_msw_segments = UDQSet::enumerateItems(context.segments());
    if (this->selector.empty()) {
        auto res = UDQSet::segments(string_value, all_msw_segments);

        auto index = std::size_t{0};
        for (const auto& ms_well : all_msw_segments) {
            for (const auto& segment : ms_well.numbers) {
                res.assign(index++, context.get_segment_var(ms_well.name,
                                                            string_value,
                                                            segment));
            }
        }

        return res;
    }

    const auto selected_segments = context.segments(this->selector);
    if (selected_segments.empty()) {
        // No matching segments.  Could be because the 'selector' only
        // applies to MS wells that are no yet online, or because the
        // segments don't yet exist.
        return UDQSet::empty(string_value);
    }
    else if (selected_segments.isScalar()) {
        // Selector matches a single segment in a single MS well.
        const auto segSet = selected_segments.segments(0);
        const auto well = std::string { segSet.well() };
        return UDQSet::scalar(string_value, context.get_segment_var(well, string_value, *segSet.begin()));
    }

    // If we get here, the selector matches at least one segment in at least
    // one MS well.
    auto res = UDQSet::segments(string_value, all_msw_segments);

    const auto numWells = selected_segments.numWells();
    for (auto wellID = 0*numWells; wellID < numWells; ++wellID) {
        const auto segSet = selected_segments.segments(wellID);
        const auto well = std::string { segSet.well() };
        for (const auto& segment : segSet) {
            res.assign(well, segment, context.get_segment_var(well, string_value, segment));
        }
    }

    return res;
}

UDQSet
UDQASTNode::eval_region_expression(const std::string& string_value,
                                   const UDQContext&  context) const
{
    const auto selected_region_sets = context.regions(string_value, this->selector);

    if (selected_region_sets.empty()) {
        // No matching region sets.  Could be because the 'selector' only
        // applies to undefined region sets, or because the regions don't
        // exist in the pertinent region set.
        return UDQSet::empty(string_value);
    }
    else if (selected_region_sets.isScalar()) {
        // Selector matches a single segment in a single MS well.
        const auto regIxRange = selected_region_sets.regions(0);
        const auto regSet = std::string { regIxRange.regionSet() };
        return UDQSet::scalar(string_value, context.get_region_var(regSet, string_value,
                                                                   *regIxRange.begin()));
    }

    // If we get here, the selector matches at least one region in at least
    // one region set.
    auto res = UDQSet::regions(string_value, UDQSet::enumerateItems(context.regions()));

    const auto numRegSets = selected_region_sets.numRegionSets();
    for (auto regSetIx = 0*numRegSets; regSetIx < numRegSets; ++regSetIx) {
        const auto regIxRange = selected_region_sets.regions(regSetIx);
        const auto regSet = std::string { regIxRange.regionSet() };
        for (const auto& regIx : regIxRange) {
            res.assign(regSet, regIx, context.get_region_var(regSet, string_value, regIx));
        }
    }

    return res;
}

UDQSet
UDQASTNode::eval_scalar_function(const UDQVarType  target_type,
                                 const UDQContext& context) const
{
    const auto& string_value = std::get<std::string>(this->value);
    const auto& udqft = context.function_table();

    const auto& func = dynamic_cast<const UDQScalarFunction&>(udqft.get(string_value));

    return func.eval(this->left->eval(target_type, context));
}

UDQSet
UDQASTNode::eval_elemental_unary_function(const UDQVarType  target_type,
                                          const UDQContext& context) const
{
    const auto& string_value = std::get<std::string>(this->value);
    const auto func_arg = this->left->eval(target_type, context);

    const auto& udqft = context.function_table();
    const auto& func = dynamic_cast<const UDQUnaryElementalFunction&>(udqft.get(string_value));

    return func.eval(func_arg);
}

UDQSet
UDQASTNode::eval_binary_function(const UDQVarType  target_type,
                                 const UDQContext& context) const
{
    const auto left_arg = this->left->eval(target_type, context);
    const auto right_arg = this->right->eval(target_type, context);
    const auto& string_value = std::get<std::string>(this->value);

    const auto& udqft = context.function_table();
    const auto& func = dynamic_cast<const UDQBinaryFunction&>(udqft.get(string_value));

    return func.eval(left_arg, right_arg);
}

UDQSet
UDQASTNode::eval_number(const UDQVarType  target_type,
                        const UDQContext& context) const
{
    const auto dummy_name = std::string { "DUMMY" };
    const auto numeric_value = std::get<double>(this->value);

    switch (target_type) {
    case UDQVarType::WELL_VAR:
        return UDQSet::wells(dummy_name, context.wells(), numeric_value);

    case UDQVarType::GROUP_VAR:
        return UDQSet::groups(dummy_name, context.nonFieldGroups(), numeric_value);

    case UDQVarType::SEGMENT_VAR:
        return UDQSet::segments(dummy_name,
                                UDQSet::enumerateItems(context.segments()),
                                numeric_value);

    case UDQVarType::SCALAR:
        return UDQSet::scalar(dummy_name, numeric_value);

    case UDQVarType::FIELD_VAR:
        return UDQSet::field(dummy_name, numeric_value);

    default:
        throw std::invalid_argument {
            "Unsupported target_type: " + std::to_string(static_cast<int>(target_type))
        };
    }
}

UDQSet
UDQASTNode::eval_table_lookup(const UDQVarType target_type,
                              const std::string& string_value,
                              const UDQContext& context) const
{
    switch (target_type) {
    case UDQVarType::FIELD_VAR:
        return eval_table_lookup_field(string_value, context);
    case UDQVarType::GROUP_VAR:
        return eval_table_lookup_group(string_value, context);
    case UDQVarType::SEGMENT_VAR:
        return eval_table_lookup_segment(string_value, context);
    case UDQVarType::WELL_VAR:
        return eval_table_lookup_well(string_value, context);
    default:
        throw std::invalid_argument {
            "Unsupported target_type: " + std::to_string(static_cast<int>(target_type))
        };
    }
}

UDQSet
UDQASTNode::eval_table_lookup_field(const std::string& string_value,
                                    const UDQContext& context) const
{
    const UDT& udt = context.get_udt(string_value);
    const auto xvar = context.get(this->selector[0]);
    double val = 0.0;
    if (xvar.has_value()) {
        val = udt(*xvar);
    }

    return UDQSet::scalar("dummy", val);
}

UDQSet
UDQASTNode::eval_table_lookup_group(const std::string& string_value,
                                    const UDQContext& context) const
{
    const UDT& udt = context.get_udt(string_value);

    const auto groups = context.nonFieldGroups();

    UDQSet result = UDQSet::groups("dummy", groups);
    for (const auto& group : groups) {
        const auto xvar = context.get_group_var(group, this->selector[0]);
        if (xvar.has_value()) {
            result.assign(group, udt(*xvar));
        }
    }

    return result;
}

UDQSet
UDQASTNode::eval_table_lookup_segment(const std::string& string_value,
                                      const UDQContext& context) const
{
    const UDT& udt = context.get_udt(string_value);

    const auto all_msw_segments = UDQSet::enumerateItems(context.segments());
    UDQSet result = UDQSet::segments("dummy", all_msw_segments);

    for (const auto& ms_well : all_msw_segments) {
        for (const auto& segment : ms_well.numbers) {
            const auto xvar = context.get_segment_var(ms_well.name,
                                                      this->selector.front(),
                                                      segment);

            if (xvar.has_value()) {
                result.assign(ms_well.name, segment, udt(*xvar));
            }
        }
    }

    return result;
}

UDQSet
UDQASTNode::eval_table_lookup_well(const std::string& string_value,
                                   const UDQContext& context) const
{
    const UDT& udt = context.get_udt(string_value);
    UDQSet result = UDQSet::wells("dummy", context.wells());
    for (const auto& well : context.wells()) {
        const auto xvar = context.get_well_var(well, this->selector[0]);
        if (xvar.has_value()) {
            result.assign(well, udt(*xvar));
        }
    }

    return result;
}

void UDQASTNode::func_tokens(std::set<UDQTokenType>& tokens) const
{
    tokens.insert(this->type);

    if (this->left) {
        this->left->func_tokens(tokens);
    }

    if (this->right) {
        this->right->func_tokens(tokens);
    }
}

void UDQASTNode::populateRequiredObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    const auto vectorType = UDQ::targetType(std::get<std::string>(this->value));

    if (this->selector.empty() && (vectorType != UDQVarType::REGION_VAR)) {
        // No specific objects identified for this node.  Will generate an
        // appropriately sized UDQ set and populate all elements of the set.
        // Nothing to do.
        //
        // Note empty selector exception in the case of REGION_VARs.  Region
        // level vectors encode the region set into the vector name, so we
        // always have to run the full requisite object analysis for region
        // level vectors.
        return;
    }

    // If we get here, then the node identifies one or more specific objects
    // for which to generate values in the resulting UDQ set.  Tell caller
    // about these objects.

    switch (vectorType) {
    case UDQVarType::REGION_VAR:
        this->populateRequiredRegionObjects(objects);
        break;

    case UDQVarType::SEGMENT_VAR:
        this->populateRequiredSegmentObjects(objects);
        break;

    case UDQVarType::WELL_VAR:
        this->populateRequiredWellObjects(objects);
        break;

    case UDQVarType::GROUP_VAR:
        this->populateRequiredGroupObjects(objects);
        break;

    default:
        // Nothing to do since either FIELD_VAR or unsupported.
        break;
    }
}

void UDQASTNode::populateRequiredGroupObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    // The selector is a group name like 'G' or a group name root like 'G*'.
    assert (this->selector.size() == 1);

    objects.groups.insert(this->selector.front());
}

namespace {
    std::optional<int> parsePositiveInt(std::string_view s)
    {
        auto result = 0;
        auto [ptr, ec] { std::from_chars(s.data(), s.data() + s.size(), result) };

        if ((ec == std::errc{}) && (ptr == s.data() + s.size()) && (result > 0)) {
            // s is "7" or similar.
            return { result };
        }

        return {};
    }
} // Anonymous namespace

void UDQASTNode::populateRequiredRegionObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    const auto split = std::string::size_type{5};

    const auto& vector = std::get<std::string>(this->value);
    if (vector.size() <= split) {
        // No specific region set, meaning this vector applies to all region
        // sets.  No requisite objects to enumerate.
        return;
    }

    // If we get here, the vector is something specific like ROPR_ABC or
    // RPR__NUM.  Split off the region set portion to form a 'FIP*' array
    // name, and insert any specific region IDs into the object list.
    //
    // This is similar to, though not quite the same as, the
    // RegionSetMatcher algorithm used during node evaluation.
    auto& regNo = objects.regions["FIP" + vector.substr(split)];

    for (const auto& regID : this->selector) {
        if (const auto id = parsePositiveInt(regID); id.has_value()) {
            regNo.insert(*id);
        }
    }
}

void UDQASTNode::populateRequiredSegmentObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    auto& segNo = objects.msWells[this->selector.front()];

    for (auto i = 0*this->selector.size() + 1; i < this->selector.size(); ++i) {
        if (const auto id = parsePositiveInt(this->selector[i]); id.has_value()) {
            segNo.insert(*id);
        }
    }
}

void UDQASTNode::populateRequiredWellObjects(UDQ::RequisiteEvaluationObjects& objects) const
{
    // The selector is a well name like 'P', a well name template like 'P*',
    // a well list like '*P', or a well list template like '*P*'.
    assert (this->selector.size() == 1);

    objects.wells.insert(this->selector.front());
}

UDQASTNode operator*(const UDQASTNode&lhs, double sign_factor)
{
    UDQASTNode prod = lhs;

    prod.scale(sign_factor);

    return prod;
}

UDQASTNode operator*(double lhs, const UDQASTNode& rhs)
{
    return rhs * lhs;
}

} // namespace Opm
