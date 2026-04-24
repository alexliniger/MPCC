#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "ADCodeGen/ad_dynamics.h"
#include "Config/config.h"
#include "Config/params.h"
#include "Config/types.h"

using json = nlohmann::json;
using mpcc::ADDynamics;
using mpcc::IntegratorType;
using mpcc::NU;
using mpcc::NX;
using mpcc::si_index;
using mpcc::StateVector;

bool generateLibraries(const std::string &configPath,
                       const std::string &modelPath,
                       const std::string &outputPath, double Ts, int n_steps) {
  std::filesystem::create_directories(outputPath);
  std::cout << "Generating libraries to: " << outputPath
            << " (Steps: " << n_steps << ")..." << std::endl;
  ADDynamics ad_dyn(Ts, modelPath, outputPath);

  std::cout << "  Generating RK4..." << std::endl;
  ad_dyn.genLibraryIntegrator(IntegratorType::RK4, n_steps);
  std::cout << "  Generating Euler..." << std::endl;
  ad_dyn.genLibraryIntegrator(IntegratorType::ForwardEuler, n_steps);

  ad_dyn.genLibraryGetF();
  ad_dyn.genLibraryTireFront();
  ad_dyn.genLibraryTireRear();

  // Generate hash of config and model files
  std::hash<std::string> hasher;
  std::ifstream iConfig(configPath);
  if (!iConfig.is_open()) {
    std::cerr << "Error: Could not open config file for hashing: " << configPath
              << std::endl;
    return false;
  }
  json jsonConfig;
  iConfig >> jsonConfig;
  size_t config_hash = hasher(jsonConfig.dump());

  std::ifstream iModel(modelPath);
  if (!iModel.is_open()) {
    std::cerr << "Error: Could not open model file: " << modelPath << std::endl;
    return false;
  }
  json jsonModel;
  iModel >> jsonModel;
  size_t model_hash = hasher(jsonModel.dump());

  std::ofstream hashFile(outputPath + "/params_hash.txt");
  if (!hashFile.is_open()) {
    std::cerr << "Error: Could not open " << outputPath
              << "/params_hash.txt for writing." << std::endl;
    return false;
  }
  hashFile << config_hash << " " << model_hash;
  hashFile.close();
  std::cout << "Libraries generated successfully." << std::endl;
  return true;
}

void verifyGeneratedLibraries(const std::string &outputPath,
                              IntegratorType type) {
  std::string model_name = (type == IntegratorType::RK4 ? "RK4" : "Euler");
  std::string lib_name = outputPath + "/cppad_cg_" + model_name + ".so";
  std::cout << "Verifying " << model_name << " library (" << lib_name << ")..."
            << std::endl;

  std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> myLib =
      std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(lib_name);
  std::unique_ptr<CppAD::cg::GenericModel<double>> model;
  model = myLib->model(model_name);

  Eigen::Matrix<double, NX + NU, 1> x;
  x.setZero();
  x(si_index.vx) = 10.0;
  std::vector<double> xv(x.data(), x.data() + x.size());
  Eigen::MatrixXd jac = Eigen::Map<Eigen::Matrix<double, NX *(NX + NU), 1>>(
      (model->Jacobian(xv)).data());
  jac.resize(NX + NU, NX);
  Eigen::MatrixXd jac2 = jac.transpose();
  std::cout << model_name << " Jacobian:\n" << jac2 << std::endl;

  StateVector fw = Eigen::Map<StateVector>((model->ForwardZero(xv)).data());
  std::cout << model_name << " ForwardZero:\n" << fw << std::endl;

  std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> cont_dyn_lib =
      std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(
          outputPath + "/cppad_cg_f_dyn.so");
  std::unique_ptr<CppAD::cg::GenericModel<double>> model_cont_dyn;
  model_cont_dyn = cont_dyn_lib->model("f_dyn");
  StateVector f =
      Eigen::Map<StateVector>((model_cont_dyn->ForwardZero(xv)).data());
  std::cout << "f_dyn ForwardZero:\n" << f << std::endl;
}

std::string resolvePath(const std::string &filename) {
  // Common locations for Params folder relative to ADCodeGen execution
  std::vector<std::string> searchPaths = {
      "../../Params/" + filename,  // from build/
      "../Params/" + filename,     // from subproject root
      "Params/" + filename         // from repo root
  };

  for (const auto &path : searchPaths) {
    std::ifstream f(path);
    if (f.good()) return path;
  }
  return "";  // Not found
}

int main(int argc, char **argv) {
  std::string configPath = "";
  std::string modelPath = "";
  std::string outputPath = "";
  bool verify = false;

  IntegratorType verification_type = IntegratorType::RK4;
  int n_steps = 1;

  // CLI Parsing
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--config" && i + 1 < argc) {
      configPath = argv[++i];
    } else if (arg == "--model" && i + 1 < argc) {
      modelPath = argv[++i];
    } else if (arg == "--output" && i + 1 < argc) {
      outputPath = argv[++i];
    } else if (arg == "--verify") {
      verify = true;
    } else if (arg == "--integrator" && i + 1 < argc) {
      std::string int_str = argv[++i];
      if (int_str == "Euler" || int_str == "ForwardEuler") {
        verification_type = IntegratorType::ForwardEuler;
      } else {
        verification_type = IntegratorType::RK4;
      }
    } else if (arg == "--steps" && i + 1 < argc) {
      n_steps = std::stoi(argv[++i]);
    } else if (arg == "--help") {
      std::cout << "Usage: ADCodeGen [options]\n"
                << "Options:\n"
                << "  --config PATH      Path to config.json\n"
                << "  --model PATH       Path to model.json\n"
                << "  --output PATH      Directory for generated libraries\n"
                << "  --integrator TYPE  RK4 or Euler\n"
                << "  --steps N          Integration steps\n"
                << "  --verify           Run verification\n"
                << "  --help             Show help\n";
      return 0;
    }
  }

  // Resolve defaults
  if (configPath.empty()) configPath = resolvePath("config.json");
  if (configPath.empty() || !std::ifstream(configPath).good()) {
    std::cerr << "Error: Could not find config.json." << std::endl;
    return 1;
  }

  std::ifstream iConfig(configPath);
  json jsonConfig;
  iConfig >> jsonConfig;

  if (modelPath.empty()) {
    if (jsonConfig.contains("model_path")) {
      std::string pathInConfig = jsonConfig["model_path"];
      if (std::ifstream(pathInConfig).good()) {
        modelPath = pathInConfig;
      } else {
        // If path starts with "Params/", strip it to use resolvePath helper
        std::string prefix = "Params/";
        if (pathInConfig.rfind(prefix, 0) == 0) {
          modelPath = resolvePath(pathInConfig.substr(prefix.length()));
        } else {
          // Try resolving as given
          modelPath = resolvePath(pathInConfig);
        }
      }
    } else {
      modelPath = resolvePath("model.json");
    }
  }

  if (modelPath.empty() || !std::ifstream(modelPath).good()) {
    std::cerr << "Error: Could not find model file: " << modelPath << std::endl;
    return 1;
  }

  if (outputPath.empty()) {
    // Determine project root from configPath to handle relative output paths
    std::string projectRoot = "./";
    if (configPath.rfind("../../Params/", 0) == 0) {
      projectRoot = "../../";
    } else if (configPath.rfind("../Params/", 0) == 0) {
      projectRoot = "../";
    }

    if (jsonConfig.contains("adcodegen_path")) {
      std::string pathInConfig = jsonConfig["adcodegen_path"];
      // If path is relative, prepend project root
      if (pathInConfig.length() > 0 && pathInConfig[0] != '/' &&
          pathInConfig[0] != '.') {
        outputPath = projectRoot + pathInConfig;
      } else {
        outputPath = pathInConfig;
      }
    } else {
      outputPath = projectRoot + "Generated";
    }
  }

  // Load defaults from config if not overridden by CLI
  if (jsonConfig.contains("integrator") && argc == 1) {
    std::string int_str = jsonConfig["integrator"];
    if (int_str == "Euler" || int_str == "ForwardEuler") {
      verification_type = IntegratorType::ForwardEuler;
    }
  }
  if (jsonConfig.contains("n_steps") && argc == 1) {
    n_steps = jsonConfig["n_steps"];
  }

  if (!generateLibraries(configPath, modelPath, outputPath, jsonConfig["Ts"],
                         n_steps)) {
    return 1;
  }

  if (verify) {
    verifyGeneratedLibraries(outputPath, verification_type);
  }

  return 0;
}
