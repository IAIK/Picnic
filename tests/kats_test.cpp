/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <picnic.h>

#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

struct test_vector {
  std::vector<uint8_t> message;
  std::vector<uint8_t> signature;
  std::vector<uint8_t> pk;
  std::vector<uint8_t> sk;
  size_t message_length;
  size_t signature_length;
};

namespace {
  bool starts_with(const std::string& str, const std::string& prefix) {
    if (str.size() < prefix.size()) {
      return false;
    }

    return std::equal(prefix.begin(), prefix.end(), str.begin());
  }

  std::vector<uint8_t> read_hex(const std::string& data) {
    std::vector<uint8_t> res;
    res.reserve(data.size() / 2);

    for (std::size_t i = 0; i < data.size(); i += 2) {
      std::istringstream iss{data.substr(i, 2)};
      unsigned int c;
      iss >> std::hex >> c;
      res.emplace_back(c & 0xff);
    }
    return res;
  }

  std::istream& operator>>(std::istream& in, test_vector& tv) {
    std::string line;
    bool expect_data = false;

    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') {
        if (expect_data) {
          // std::cout << "Expected data." << std::endl;
          in.setstate(std::ios::failbit);
          break;
        }
        // skip empty lines and comments
        continue;
      }

      if (starts_with(line, "count = ")) {
        // skip count
        expect_data = true;
        continue;
      } else if (starts_with(line, "seed = ")) {
        // skip seed
        continue;
      } else if (starts_with(line, "mlen = ")) {
        // read message length
        std::istringstream iss{line.substr(7)};
        iss >> tv.message_length;
        if (!iss) {
          // std::cout << "Unable to parse message length." << std::endl;
          break;
        }
      } else if (starts_with(line, "smlen = ")) {
        // read signature length
        std::istringstream iss{line.substr(8)};
        iss >> tv.signature_length;
        if (!iss) {
          // std::cout << "Unable to parse signature length." << std::endl;
          break;
        }
      } else if (starts_with(line, "msg = ")) {
        // read message
        tv.message = read_hex(line.substr(6));
      } else if (starts_with(line, "pk = ")) {
        // read pk
        tv.pk = read_hex(line.substr(5));
      } else if (starts_with(line, "sk = ")) {
        // read sk
        tv.sk = read_hex(line.substr(5));
      } else if (starts_with(line, "sm = ")) {
        // read signature
        tv.signature = read_hex(line.substr(5));
        expect_data  = false;
        break;
      } else {
        // std::cout << "Do not know how to handle line: " << line << std::endl;
        in.setstate(std::ios::failbit);
        break;
      }
    }

    return in;
  }

  bool run_picnic_test(const test_vector& tv) {
    picnic_privatekey_t private_key;
    picnic_publickey_t public_key;

    int ret = picnic_read_private_key(&private_key, tv.sk.data(), tv.sk.size());
    if (ret != 0) {
      std::cout << "Unable to read private key." << std::endl;
      return false;
    }

    ret = picnic_read_public_key(&public_key, tv.pk.data(), tv.pk.size());
    if (ret != 0) {
      std::cout << "Unable to read public key." << std::endl;
      return false;
    }

    ret = picnic_validate_keypair(&private_key, &public_key);
    if (ret != 0) {
      std::cout << "Key pair does not validate." << std::endl;
      return false;
    }

    // Test vectors generated for NIST have message length and the message at the beginning.
    const size_t offset = tv.message.size() + sizeof(uint32_t);
    std::vector<uint8_t> tvsig{tv.signature.begin() + offset, tv.signature.end()};

    std::vector<uint8_t> signature;
    signature.resize(tvsig.size() * 2);

    // Recreate the signature
    size_t signature_len = signature.size();
    ret = picnic_sign(&private_key, tv.message.data(), tv.message.size(), signature.data(),
                      &signature_len);
    if (ret != 0) {
      std::cout << "Unable to sign." << std::endl;
      return false;
    }
    signature.resize(signature_len);

    // Verify the provided signature
    ret = picnic_verify(&public_key, tv.message.data(), tv.message.size(), tvsig.data(),
                        tvsig.size());
    if (ret != 0) {
      std::cout << "Signature does not verify." << std::endl;
      return false;
    }

    // Check if computed signature matches
    if (signature.size() != tvsig.size()) {
      std::cout << "Signature length does not match." << std::endl;
      return false;
    }
    if (signature != tvsig) {
      std::cout << "Signature does not match." << std::endl;
      return false;
    }

    return true;
  }

  bool run_test_vectors_from_file(const char* path, size_t pks, size_t sks) {
    std::ifstream in{path};
    if (!in) {
      return false;
    }

    size_t vectors_run       = 0;
    size_t vectors_succeeded = 0;
    test_vector tv;
    while (in >> tv) {
      if (tv.sk.size() != sks) {
        std::cout << "Invalid secret key length." << std::endl;
        continue;
      }
      if (tv.pk.size() != pks) {
        std::cout << "Invalid public key length." << std::endl;
        continue;
      }
      if (tv.message.size() != tv.message_length || tv.signature.size() != tv.signature_length) {
        std::cout << "Invalid message or signature length." << std::endl;
        continue;
      }

      ++vectors_run;
      vectors_succeeded += run_picnic_test(tv) ? 1 : 0;

      tv = test_vector{};
    };

    return vectors_run && vectors_succeeded == vectors_run;
  }

  int run_test(picnic_params_t param) {
    static constexpr const char* tests[] = {
        NULL,
        KATDIR "/kat_l1_fs.txt",
        KATDIR "/kat_l1_ur.txt",
        KATDIR "/kat_l3_fs.txt",
        KATDIR "/kat_l3_ur.txt",
        KATDIR "/kat_l5_fs.txt",
        KATDIR "/kat_l5_ur.txt",
        KATDIR "/kat_picnic3_l1.txt",
        KATDIR "/kat_picnic3_l3.txt",
        KATDIR "/kat_picnic3_l5.txt",
        KATDIR "/kat_l1_full.txt",
        KATDIR "/kat_l3_full.txt",
        KATDIR "/kat_l5_full.txt",
    };

    return run_test_vectors_from_file(tests[param], picnic_get_public_key_size(param),
                                      picnic_get_private_key_size(param));
  }
} // namespace

int main(int argc, char** argv) {
  if (argc == 2) {
    const picnic_params_t param = argument_to_params(argv[1]);
    if (param == PARAMETER_SET_INVALID) {
      std::cout << "ERR: invalid test idx" << std::endl;
      return 1;
    }

    const int t = run_test(param);
    if (!t) {
      std::cout << "ERR: Picnic KAT test " << picnic_get_param_name(param) << " FAILED (" << t
                << ")" << std::endl;
      return -1;
    }
    return 0;
  }

  int ret = 0;
  for (const auto s : all_supported_parameters()) {
    const int t = run_test(s);
    if (!t) {
      std::cout << "ERR: Picnic KAT test: " << picnic_get_param_name(s) << " FAILED (" << t << ")"
                << std::endl;
      ret = -1;
    }
  }

  return ret;
}
