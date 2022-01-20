#include "utils.h"

#include <algorithm>
#include <string>
#include <sstream>

namespace {
  constexpr picnic_params_t all_parameters[] = {
      Picnic_L1_FS, Picnic_L1_UR, Picnic_L3_FS, Picnic_L3_UR,   Picnic_L5_FS,   Picnic_L5_UR,
      Picnic3_L1,   Picnic3_L3,   Picnic3_L5,   Picnic_L1_full, Picnic_L3_full, Picnic_L5_full,
  };

  bool isequal(const std::string& lhs, const std::string& rhs) {
    return lhs.size() == rhs.size() &&
           std::equal(lhs.begin(), lhs.end(), rhs.begin(),
                      [](const std::string::value_type l, const std::string::value_type r) {
                        return std::tolower(l) == std::tolower(r);
                      });
  }
} // namespace

std::vector<picnic_params_t> all_supported_parameters() {
  std::vector<picnic_params_t> res;
  for (const auto param : all_parameters) {
    if (picnic_signature_size(param)) {
      res.push_back(param);
    }
  }
  return res;
}

picnic_params_t argument_to_params(const char* arg) {
  const std::string sarg{arg};

  for (const auto param : all_parameters) {
    const std::string name{picnic_get_param_name(param)};
    if (isequal(sarg, name)) {
      return param;
    }
  }

  std::istringstream iss{sarg};
  unsigned int idx;
  iss >> idx;
  if (!iss || idx == PARAMETER_SET_INVALID || idx >= PARAMETER_SET_MAX_INDEX) {
    return PARAMETER_SET_INVALID;
  }

  return static_cast<picnic_params_t>(idx);
}