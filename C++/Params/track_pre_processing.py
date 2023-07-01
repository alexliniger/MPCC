import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev


def project_to_spline(pos, tck, s_start=0, comp_initial_guess=True):
    if comp_initial_guess:
        s_full = np.linspace(0, 1, 1000)
        x_full = splev(s_full, tck)
        x_full = np.array(x_full)
        dist = np.sum(
            (x_full - np.repeat(pos[:, np.newaxis], 1000, axis=1)) ** 2, axis=0
        )
        i_min = np.argmin(dist)
        s_opt = s_full[i_min]
    else:
        s_opt = s_start
    error = None
    for i in range(20):
        pos_spline = splev(s_opt, tck)
        ds_spline = splev(s_opt, tck, der=1)
        dds_spline = splev(s_opt, tck, der=2)
        diff = pos_spline - pos
        jac = 2.0 * diff[0] * ds_spline[0] + 2.0 * diff[1] * ds_spline[1]
        hessian = (
            2.0 * ds_spline[0] ** 2
            + 2.0 * diff[0] * dds_spline[0]
            + 2.0 * ds_spline[1] ** 2
            + 2.0 * diff[1] * dds_spline[1]
        )

        s_opt -= jac / hessian
        s_opt = s_opt % 1
        error = np.linalg.norm(diff)

    return s_opt, error


def interpolate_track(path, json_file):
    # Load data from JSON file
    with open(os.path.join(path, json_file), "r") as file:
        track = json.load(file)

    # Perform 2D interpolation for center/race line
    x = np.array(track["X"])
    y = np.array(track["Y"])
    # It is often necessary to slightly smooth tracks such that the curvature
    # is not too noisy. For the smoothing we use the weight w of each point
    # and set it to the inverse of the standard deviation, since this is not
    # really known "w becomes a tuning variable". Following the guide of splprep
    # this allows to set the smoothness weight to the number of elements that
    # are interpolated.
    track_uncertainty = 0.025  # standard deviation of track points in m
    w = 1 / track_uncertainty * np.ones(len(x))
    tck, u = splprep([x, y], s=len(x), w=w)

    # Perform 2D interpolation for outer line
    x_outer = np.array(track["X_o"])
    y_outer = np.array(track["Y_o"])
    w = 1 / track_uncertainty * np.ones(len(x_outer))
    tck_outer, u_outer = splprep([x_outer, y_outer], s=len(x_outer), w=w)

    # Perform 2D interpolation for outer line
    x_inner = np.array(track["X_i"])
    y_inner = np.array(track["Y_i"])
    w = 1 / track_uncertainty * np.ones(len(x_inner))
    tck_inner, u_innerr = splprep([x_inner, y_inner], s=len(x_inner), w=w)

    # Generate points for interpolation
    n_elements = 1000
    s = np.linspace(0, 1, n_elements)

    # Evaluate interpolated points
    x_center_inter = splev(s, tck)
    x_outer_inter = splev(s, tck_outer)
    x_inner_inter = splev(s, tck_inner)

    # compute distance to inner and outer border
    # this assumes that the "center" line is never outside the borders
    n_inner = []
    n_outer = []
    for i in range(len(s)):
        s_i_inner, n_i_inner = project_to_spline(
            np.array(x_center_inter)[:, i], tck_inner
        )
        s_i_outer, n_i_outer = project_to_spline(
            np.array(x_center_inter)[:, i], tck_outer
        )
        n_inner.append(n_i_inner)
        n_outer.append(n_i_outer)
    n_outer = np.array(n_outer)
    n_inner = np.array(n_inner)
    # compute curvature
    ds_center = splev(s, tck, der=1)
    dds_center = splev(s, tck, der=2)
    curvature = (ds_center[0] * dds_center[1] - ds_center[1] * dds_center[0]) / (
        ds_center[0] ** 2 + ds_center[1] ** 2
    ) ** 1.5

    # compute metric track length
    s_long = np.linspace(0, 1, 10000)
    x_long = splev(s_long, tck)
    x_long_diff = np.diff(x_long)
    track_length = np.sum(np.linalg.norm(x_long_diff, axis=0))

    s_metric = s * track_length

    # compute safe speed
    # this uses the idea of Safe Motion Planning for Autonomous Driving
    # using an Adversarial Road Model, there are 3 variables:
    # maximum lateral acceleration
    a_max = 15.0
    # maximum longitudinal velocity
    v_max = 250 / 3.6
    # lookahead to be considered, depends on the tack layout
    # rule of thumb, slightly more than distance to slow down from
    # top speed to min speed on track
    lookahead = 200  # lookahead in m
    lookahead_points = int(lookahead / (track_length / n_elements))
    safe_speed = []
    for i in range(n_elements):
        indices = np.arange(i, i + lookahead_points) % n_elements
        curvature_max = np.max(np.abs(curvature[indices]))
        safe_speed.append(np.minimum(np.sqrt(a_max / curvature_max), v_max))
    safe_speed = np.array(safe_speed)

    track_processed = {
        "X": x_center_inter[0].tolist(),
        "Y": x_center_inter[1].tolist(),
        "X_o": x_outer_inter[0].tolist(),
        "Y_o": x_outer_inter[1].tolist(),
        "X_i": x_inner_inter[0].tolist(),
        "Y_i": x_inner_inter[1].tolist(),
        "s": s_metric.tolist(),
        "n_left": n_outer.tolist(),
        "n_right": (-n_inner).tolist(),
        "velocity": safe_speed.tolist(),
    }

    with open(os.path.join(path, "track_processed.json"), "w") as file:
        json.dump(track_processed, file)

    # Create a 2D plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Plot original points
    ax.plot(x, y, "-r", label="Original Points")
    ax.plot(x_center_inter[0], x_center_inter[1], "-b", label="Interpolated Points")

    ax.plot(x_inner, y_inner, "-r")
    ax.plot(x_inner_inter[0], x_inner_inter[1], "-b")

    ax.plot(x_outer, y_outer, "-r")
    ax.plot(x_outer_inter[0], x_outer_inter[1], "-b")

    ax.axis("equal")
    ax.legend()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(s_metric, -n_inner)
    ax.plot(s_metric, n_outer)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(s_metric, curvature)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(s_metric, safe_speed)

    # Display the plot
    plt.show()


# Call the function with the JSON file path
interpolate_track(path="C++/Params/", json_file="track.json")
