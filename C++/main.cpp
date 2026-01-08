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
#include <string>

#include <nlohmann/json.hpp>

#include "Config/track.h"
#include "MPC/mpc.h"
#include "Model/integrator.h"
#include "Plotting/plotting.h"
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
  mpc.setTrack(track_full);
  const double phi_0 =
      std::atan2(track_xy.Y(1) - track_xy.Y(0), track_xy.X(1) - track_xy.X(0));
  mpcc::State x0 = {
      track_xy.X(0), track_xy.Y(0), phi_0, jsonConfig["v0"], 0.0, 0.0, 0.0,
      0.0,           0.0,           0.0,   jsonConfig["v0"]};
  for (int i = 0; i < jsonConfig["n_sim"]; i++) {
    mpcc::MPCReturn mpc_sol = mpc.runMPC(x0);
    // Use the MPC prediction as sim step
    // x0 = mpc_sol.mpc_horizon[1].xk;
    // Use ODE integrator
    x0 = integrator.simTimeStep(x0, mpc_sol.u0, jsonConfig["Ts"]);
    x0.unwrap(track_full.s(track_full.s.size() - 1));
    log.push_back(mpc_sol);
  }

  double mean_time = 0.0;
  double max_time = 0.0;
  for (const mpcc::MPCReturn &log_i : log) {
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
