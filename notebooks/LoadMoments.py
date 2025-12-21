import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd

def read_solution_file(filename):
    """
    Reads the Trixi-like solution file into a dict:
    {
        timestep: {
            't': float,
            'x': np.ndarray,
            'w': np.ndarray  # shape (n_points, n_vars)
        },
        ...
    }
    """
    data = {}
    timestep = None
    x_vals = []
    u_vals = []

    with open(filename, "r") as f:
        counter = 0
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("# timestep"):
                # Save previous timestep if available
                if timestep is not None:
                    data[counter] = {
                        "t": time,
                        "x": np.array(x_vals),
                        "u": np.array(u_vals)
                    }
                    counter += 1
                # Parse new header
                parts = line.replace(",", "").split()
                timestep = int(parts[1].split("=")[-1]) if "timestep=" in line else int(parts[1])
                points = int(parts[2].split("=")[-1])
                time = float(parts[3].split("=")[-1])
                x_vals, u_vals = [], []

            elif not line.startswith("#"):
                vals = list(map(float, line.split()))
                x_vals.append(vals[0])
                u_vals.append(vals[1:])  # w^(0..n)

        # Save last timestep
        if timestep is not None:
            data[counter] = {
                "t": time,
                "x": np.array(x_vals),
                "u": np.array(u_vals)
            }
            counter += 1

    return data

def interpolate_solution(data, time_query, times, timesteps):
    """Linear interpolation in time for all w's"""
    # Clamp outside range
    if time_query <= times.min():
        return data[timesteps[0]]
    if time_query >= times.max():
        return data[timesteps[-1]]

    # indices around query
    idx_right = np.searchsorted(times, time_query)
    idx_left = idx_right - 1
    t0, t1 = times[idx_left], times[idx_right]
    u0, u1 = data[timesteps[idx_left]]["u"], data[timesteps[idx_right]]["u"]
    x0 = data[timesteps[idx_left]]["x"]

    # linear interpolation in time
    alpha = (time_query - t0) / (t1 - t0)
    w_interp = (1 - alpha) * u0 + alpha * u1

    return {"t": time_query, "x": x0, "w": w_interp}

def plot_solution_time(data, time_query, times, timesteps):
    sol = interpolate_solution(data, time_query, times, timesteps)
    x, u, t = sol["x"], sol["u"], sol["t"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

    axes[0].plot(x, u[:,0], label="$u^{(0)}$", color='blue')
    axes[0].plot(x, u[:,1], label="$u^{(1)}$", color='orange')
    axes[0].plot(x, u[:,2], label="$v^{(2)}$", color='green')
    axes[0].set_title(f"t={t:.4f}")
    axes[0].set_xlabel("X")
    axes[0].grid()
    axes[0].legend()

    axes[1].plot(x, u[:,3], label="$u^{(3)}$", color='purple')
    axes[1].plot(x, u[:,4], label="$u^{(4)}$", color='pink')
    # axes[1].plot(x, w[:,5], label="$w^{(5)}$", color='brown')
    # axes[1].plot(x, w[:,6], label="$w^{(6)}$", color='gray')
    axes[1].set_title("Higher moments")
    axes[1].set_xlabel("X")
    axes[1].grid()
    axes[1].legend()

    plt.tight_layout()
    plt.show()




def read_solution_file_bgk(filename):
    """
    Reads the Trixi-like solution file into a dict:
    {
        timestep: {
            't': float,
            'x': np.ndarray,
            'w': np.ndarray  # shape (n_points, n_vars)
        },
        ...
    }
    """
    data = {}
    timestep = None
    x_vals = []
    f_vals = []

    with open(filename, "r") as f:
        counter = 0
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("# timestep"):
                # Save previous timestep if available
                if timestep is not None:
                    data[counter] = {
                        "t": time,
                        "x": np.array(x_vals),
                        "f": np.array(f_vals)
                    }
                    counter += 1
                # Parse new header
                parts = line.replace(",", "").split()
                timestep = int(parts[1].split("=")[-1]) if "timestep=" in line else int(parts[1])
                points = int(parts[2].split("=")[-1])
                time = float(parts[3].split("=")[-1])
                x_vals, f_vals = [], []

            elif not line.startswith("#"):
                vals = list(map(float, line.split()))
                x_vals.append(vals[0])
                f_vals.append(vals[1:])  # w^(0..n)

        # Save last timestep
        if timestep is not None:
            data[counter] = {
                "t": time,
                "x": np.array(x_vals),
                "f": np.array(f_vals)
            }
            counter += 1

    return data

def interpolate_solution_bgk(data, time_query, times, timesteps):
    """Linear interpolation in time for all w's"""
    # Clamp outside range
    if time_query <= times.min():
        return data[timesteps[0]]
    if time_query >= times.max():
        return data[timesteps[-1]]

    # indices around query
    idx_right = np.searchsorted(times, time_query)
    idx_left = idx_right - 1
    t0, t1 = times[idx_left], times[idx_right]
    f0, f1 = data[timesteps[idx_left]]["f"], data[timesteps[idx_right]]["f"]
    x0 = data[timesteps[idx_left]]["x"]

    # linear interpolation in time
    alpha = (time_query - t0) / (t1 - t0)
    f_interp = (1 - alpha) * f0 + alpha * f1

    return {"t": time_query, "x": x0, "f": f_interp}