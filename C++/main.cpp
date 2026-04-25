// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <string>

#include <nlohmann/json.hpp>

#include "Config/track.h"
#include "MPC/mpc.h"
#include "Model/integrator.h"
#include "Model/sim_model.h"
#include "Plotting/plotting.h"
#include "Spline/arc_length_spline.h"
#include "Spline/boost_splines.h"
using json = nlohmann::json;

int main() {
  std::string go_to_path = "../";
  std::ifstream iConfig(go_to_path + "Params/config.json");
  json jsonConfig;
  iConfig >> jsonConfig;

  mpcc::PathToJson json_paths{
      go_to_path + std::string(jsonConfig["model_path"]),
      go_to_path + std::string(jsonConfig["cost_path"]),
      go_to_path + std::string(jsonConfig["bounds_path"]),
      go_to_path + std::string(jsonConfig["track_path"]),
      go_to_path + std::string(jsonConfig["normalization_path"]),
      go_to_path + std::string(jsonConfig["adcodegen_path"])};

  mpcc::Track track = mpcc::Track(json_paths.track_path);
  mpcc::TrackPos track_xy = track.getTrack();
  mpcc::TrackFull track_full = track.getTrackFull();

  mpcc::Integrator integrator = mpcc::Integrator(jsonConfig["Ts"], json_paths);
  mpcc::Plotting plotter =
      mpcc::Plotting(jsonConfig["Ts"], json_paths, track_full);

  std::list<mpcc::MPCReturn> log;
  mpcc::MPC mpc(jsonConfig["n_sqp"], jsonConfig["n_reset"],
                jsonConfig["sqp_mixing"], jsonConfig["Ts"], json_paths);
  std::cout << "MPC initialized" << std::endl;
  mpc.setTrack(track_full);
  std::cout << "Track set" << std::endl;
  const double phi_0 =
      std::atan2(track_xy.Y(1) - track_xy.Y(0), track_xy.X(1) - track_xy.X(0));
  mpcc::State x0 = {
      track_xy.X(0), track_xy.Y(0), phi_0, jsonConfig["v0"], 0.0, 0.0, 0.0,
      0.0,           0.0,           0.0,   jsonConfig["v0"]};
  std::cout << "x0 initialized" << std::endl;
  mpcc::ArcLengthSpline track_spline(json_paths);
  track_spline.gen2DSpline(track_xy.X, track_xy.Y);
  std::unique_ptr<mpcc::ISimModel> sim_plant;
  if (jsonConfig.contains("sim_model_type") &&
      jsonConfig["sim_model_type"] == "four_wheel") {
    sim_plant = std::make_unique<mpcc::FourWheelSimModel>(
        jsonConfig["Ts"], json_paths.param_path);
    std::cout << "FourWheelSimModel initialized" << std::endl;
  } else {
    sim_plant = std::make_unique<mpcc::BicycleSimModel>(jsonConfig["Ts"],
                                                        json_paths.param_path);
    std::cout << "BicycleSimModel initialized" << std::endl;
  }

  // Initialize SimState from initial State
  sim_plant->initialize(x0);
  std::cout << "SimState initialized" << std::endl;

  for (int i = 0; i < jsonConfig["n_sim"]; i++) {
    std::cout << "Running MPC step " << i << std::endl;
    mpcc::MPCReturn mpc_sol = mpc.runMPC(x0);
    if (jsonConfig["use_sim_model"]) {
      // Sim using seperate model
      sim_plant->step(mpc_sol.u0, jsonConfig["Ts"]);
      // Project the x_sim position
      mpcc::State x_temp_proj = x0;
      x_temp_proj.X = sim_plant->getX();
      x_temp_proj.Y = sim_plant->getY();
      double s_proj = track_spline.projectOnSpline(x_temp_proj);

      double vs_new =
          x0.vs + mpc_sol.u0.dVs * static_cast<double>(jsonConfig["Ts"]);
      x0 = sim_plant->getMPCCState(s_proj, vs_new);
    } else {
      // Use the MPC prediction as sim step
      x0 = mpc_sol.mpc_horizon[1].xk;
    }

    x0.unwrap(track_full.s(track_full.s.size() - 1));
    log.push_back(mpc_sol);
  }

  double mean_time = 0.0;
  double max_time = 0.0;
  for (const mpcc::MPCReturn& log_i : log) {
    mean_time += log_i.time_total;
    if (log_i.time_total > max_time) max_time = log_i.time_total;
  }
  std::cout << "mean nmpc time "
            << mean_time / static_cast<double>(jsonConfig["n_sim"])
            << std::endl;
  std::cout << "max nmpc time " << max_time << std::endl;

  plotter.plotRun(log, track_xy);

  return 0;
}
