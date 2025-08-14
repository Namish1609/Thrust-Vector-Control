#!/usr/bin/env python3
"""
rocket_sim_with_czml.py

Runs the rocket simulation (PID in degrees) with initial tilt = 0 and exports
a Cesium CZML (trajectory_kennedy.czml) centered at Kennedy Space Center.

Usage:
    python3 rocket_sim_with_czml.py
"""

import json
from datetime import datetime, timedelta
import math
import numpy as np
from math import sin, cos, radians, degrees, atan2, sqrt
import matplotlib.pyplot as plt  # pyright: ignore[reportMissingModuleSource]
from scipy.integrate import solve_ivp  # pyright: ignore[reportMissingImports]
import subprocess
import os
import http.server
import socketserver
import webbrowser
import json
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# -------------------------
# PARAMETERS (unchanged)
# -------------------------
params = {
    "mass_initial": 0.9, "prop_mass": 0.28, "length": 0.90, "diameter": 0.06,
    "Cd": 0.75, "A": np.pi*(0.06/2)**2, "rho": 1.225, "Cm_alpha": -0.03, "S_ref": 0.02,
    "cp_base": 0.02, "cp_shift_burn": -0.004, "cp_shift_coast": 0.010,
    "burn_time": 1.5, "thrust_peak": 40.0,
    "servo_speed_s_per_60deg": 0.12, "servo_max_angle_deg": 7.0,
    "rw_max_torque": 1.5, "rw_inertia": 0.002, "rw_motor_tau": 0.02,
    "mass_imbalance_g": 10.0, "cg_offset_radial": 0.02,
    # use firmware PID gains (tuned for degrees)
    "pid_pitch": {"kp": 1.15, "ki": 0.025, "kd": 0.045},
    "pid_yaw":   {"kp": 1.15, "ki": 0.025, "kd": 0.045},
    "pid_roll":  {"kp": 2.2, "ki": 0.012, "kd": 0.03},
    "control_dt": 0.003, "sensor_delay": 0.004,
    "wind_gust": {"time": 0.9, "duration": 0.12, "force_x": 01.0, "force_y": 0.60},
    "initial_tilt_pitch_deg": 45.0, "initial_tilt_yaw_deg": 0.0, "initial_roll_deg": 0.0,
    "t_final": 12.0, "dt": 0.002
}

# Derived
m0 = params["mass_initial"]; mp = params["prop_mass"]; burn_time = params["burn_time"]
length = params["length"]; r = params["diameter"]/2; A = params["A"]
I_long = 0.5 * m0 * r*r; I_trans = (1/12.0) * m0 * (3*r*r + length*length)
servo_tau = params["servo_speed_s_per_60deg"] / 2.2
servo_limit = params["servo_max_angle_deg"]

# -------------------------
# Helper classes
# -------------------------
class Kalman1D:
    def __init__(self, q=0.001, r=0.03, x0=0.0, p0=1.0):
        self.q = q; self.r = r; self.x = x0; self.p = p0
    def update(self, z, dt):
        self.p += self.q * dt
        k = self.p / (self.p + self.r)
        self.x += k * (z - self.x)
        self.p *= (1 - k)
        return self.x

class PID:
    def __init__(self,kp,ki,kd, out_min=-999, out_max=999):
        self.kp=kp; self.ki=ki; self.kd=kd
        self.i=0.0; self.last=None
        self.out_min=out_min; self.out_max=out_max
    def update(self,setpoint,meas,dt):
        err = setpoint - meas
        self.i += err * dt
        d = 0.0 if (self.last is None) else (err - self.last) / dt
        self.last = err
        out = self.kp*err + self.ki*self.i + self.kd*d
        out = max(self.out_min, min(self.out_max, out))
        if out == self.out_min or out == self.out_max:
            self.i *= 0.9
        return out
    def reset(self):
        self.i=0.0; self.last=None

# instantiate controllers + filters
pid_pitch = PID(**params["pid_pitch"])
pid_yaw   = PID(**params["pid_yaw"])
pid_roll  = PID(**params["pid_roll"])

kal_pitch = Kalman1D(q=0.001, r=0.03, x0=math.radians(params["initial_tilt_pitch_deg"]))
kal_yaw   = Kalman1D(q=0.001, r=0.03, x0=math.radians(params["initial_tilt_yaw_deg"]))
kal_alt   = Kalman1D(q=0.1, r=1.0, x0=0.0)

# -------------------------
# Kinematics / dynamics
# -------------------------
def body_to_world(vec, roll, pitch, yaw):
    cr = cos(roll); sr = sin(roll)
    cp = cos(pitch); sp = sin(pitch)
    cy = cos(yaw); sy = sin(yaw)
    R = np.array([
        [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr],
        [sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr],
        [-sp,   cp*sr,             cp*cr]
    ])
    return R.dot(vec)

def thrust_time(t):
    if t <= 0: return 0.0
    if t >= burn_time: return 0.0
    return params["thrust_peak"] * np.sin(np.pi * t / burn_time)

def mass_time(t):
    if t <= 0: return m0
    if t >= burn_time: return m0 - mp
    return m0 - (mp * (t / burn_time))

def cp_offset_time(t):
    if t <= 0: return params["cp_base"]
    if t <= burn_time:
        frac = t / burn_time
        return params["cp_base"] + params["cp_shift_burn"] * (1 - frac)
    else:
        dt = t - burn_time; tau = 1.0
        shift = params["cp_shift_coast"] * (1 - np.exp(-dt / tau))
        return params["cp_base"] + shift

def dynamics(t, state, controller_state):
    x, y, z, vx, vy, vz, pitch, pitch_rate, yaw, yaw_rate, roll, roll_rate, s_pitch, s_yaw, rw_cmd_torque, rw_torque = state
    m = mass_time(t); T = thrust_time(t); g = 9.80665
    v = sqrt(vx*vx + vy*vy + vz*vz)
    drag = 0.5 * params["rho"] * params["Cd"] * A * v*v if v>1e-6 else 0.0
    drag_x = -drag * (vx/(v+1e-9)); drag_y = -drag * (vy/(v+1e-9)); drag_z = -drag * (vz/(v+1e-9))
    gp = np.radians(s_pitch); gy = np.radians(s_yaw)
    # thrust along body +Z with small gimbals: main axis is +Z
    tb_x = T * sin(gy)
    tb_y = -T * sin(gp)
    tb_z = T * cos(gp) * cos(gy)
    thr_body = np.array([tb_x, tb_y, tb_z])
    thr_world = body_to_world(thr_body, roll, pitch, yaw)
    Tx, Ty, Tz = thr_world[0], thr_world[1], thr_world[2]
    gw = params["wind_gust"]
    gust_x = gw["force_x"] if (t>=gw["time"] and t<=gw["time"]+gw["duration"]) else 0.0
    gust_y = gw["force_y"] if (t>=gw["time"] and t<=gw["time"]+gw["duration"]) else 0.0
    ax = (Tx + drag_x + gust_x) / m
    ay = (Ty + drag_y + gust_y) / m
    az = (Tz + drag_z - m*g) / m
    alpha_pitch = pitch - (0.0 if v<1e-6 else atan2(vx, vz))
    alpha_yaw   = yaw   - (0.0 if v<1e-6 else atan2(vy, vz))
    q_dyn = 0.5 * params["rho"] * v*v
    cp = cp_offset_time(t)
    M_pitch = q_dyn * params["S_ref"] * params["Cm_alpha"] * alpha_pitch * cp
    M_yaw   = q_dyn * params["S_ref"] * params["Cm_alpha"] * alpha_yaw * cp
    torque_from_tvc_pitch = T * sin(gp) * cp
    torque_from_tvc_yaw   = T * sin(gy) * cp
    torque_imb_pitch = 0.0; torque_imb_yaw = 0.0; torque_imb_roll = 0.0
    dpitch_rate = (M_pitch + torque_from_tvc_pitch + torque_imb_pitch) / I_trans
    dyaw_rate   = (M_yaw   + torque_from_tvc_yaw   + torque_imb_yaw) / I_trans
    droll_rate  = ( -rw_torque + torque_imb_roll ) / I_long
    s_cmd_pitch = controller_state.get("s_cmd_pitch", 0.0)
    s_cmd_yaw   = controller_state.get("s_cmd_yaw", 0.0)
    ds_pitch = (s_cmd_pitch - s_pitch) / servo_tau
    ds_yaw   = (s_cmd_yaw - s_yaw) / servo_tau
    rw_cmd = controller_state.get("rw_cmd", 0.0)
    rw_cmd = max(-params["rw_max_torque"], min(params["rw_max_torque"], rw_cmd))
    drw_torque = (rw_cmd - rw_torque) / params["rw_motor_tau"]
    dx = vx; dy = vy; dz = vz; dvx = ax; dvy = ay; dvz = az
    dpitch = pitch_rate; dyaw = yaw_rate; droll = roll_rate
    return [dx, dy, dz, dvx, dvy, dvz, dpitch, dpitch_rate, dyaw, dyaw_rate, droll, droll_rate, ds_pitch, ds_yaw, 0.0, drw_torque]

def controller_step(t, state, controller_state, dt):
    x, y, z, vx, vy, vz, pitch, pitch_rate, yaw, yaw_rate, roll, roll_rate, s_pitch, s_yaw, rw_cmd_torque, rw_torque = state
    # measurements (angles in radians internally)
    meas_pitch = pitch + np.random.normal(0.0, np.radians(0.05))
    meas_yaw   = yaw   + np.random.normal(0.0, np.radians(0.05))
    meas_alt   = z     + np.random.normal(0.0, 0.1)
    # filter (still in radians)
    filt_pitch_rad = kal_pitch.update(meas_pitch, dt)
    filt_yaw_rad   = kal_yaw.update(meas_yaw, dt)
    filt_alt       = kal_alt.update(meas_alt, dt)
    # convert filtered radians -> degrees to feed PID (matches firmware units)
    filt_pitch_deg = degrees(filt_pitch_rad)
    filt_yaw_deg   = degrees(filt_yaw_rad)
    # roll rate to deg/s
    meas_roll_rate_dps = degrees(roll_rate)
    # PID (degrees)
    pitch_pid_out = pid_pitch.update(0.0, filt_pitch_deg, dt)
    yaw_pid_out   = pid_yaw.update(0.0, filt_yaw_deg, dt)
    # servo commands limited in degrees
    s_cmd_pitch = np.clip(pitch_pid_out, -servo_limit, servo_limit)
    s_cmd_yaw   = np.clip(yaw_pid_out,   -servo_limit, servo_limit)
    # roll PID (rate controller): input deg/s, output torque (mapped)
    rw_out = pid_roll.update(0.0, meas_roll_rate_dps, dt)
    rw_out = np.clip(rw_out, -params["rw_max_torque"], params["rw_max_torque"])
    controller_state["s_cmd_pitch"] = s_cmd_pitch
    controller_state["s_cmd_yaw"] = s_cmd_yaw
    controller_state["rw_cmd"] = rw_out
    return filt_pitch_rad, filt_yaw_rad, filt_alt

# -------------------------
# RUN SIMULATION
# -------------------------
init_pitch = math.radians(params["initial_tilt_pitch_deg"])
init_yaw   = math.radians(params["initial_tilt_yaw_deg"])
init_roll  = math.radians(params["initial_roll_deg"])

state0 = [0.0,0.0,0.0,   # x,y,z (m)
          0.0,0.0,0.0,   # vx,vy,vz
          init_pitch,0.0, # pitch, pitch_rate
          init_yaw,0.0,   # yaw, yaw_rate
          init_roll,0.0,  # roll, roll_rate
          0.0,0.0,        # s_pitch, s_yaw (servo positions in deg)
          0.0,0.0]        # rw_cmd_torque, rw_torque

dt = params["dt"]
times = np.arange(0.0, params["t_final"] + dt, dt)
states = np.zeros((len(times), len(state0)))
states[0,:] = state0
controller_state = {"s_cmd_pitch":0.0, "s_cmd_yaw":0.0, "rw_cmd":0.0}
events = []
flight_state = "GROUND"
boost_start = None
apogee_time = None
parachute_deployed = False
max_alt = 0.0
next_control_time = 0.0

for i in range(1, len(times)):
    t = times[i-1]
    sol = solve_ivp(lambda tt, y: dynamics(tt, y, controller_state), [t, t+dt], states[i-1,:], method='RK45', max_step=dt)
    states[i,:] = sol.y[:, -1]
    if times[i] >= next_control_time:
        filt_pitch, filt_yaw, filt_alt = controller_step(times[i] + params["sensor_delay"], states[i,:], controller_state, params["control_dt"])
        next_control_time += params["control_dt"]
    else:
        filt_pitch, filt_yaw, filt_alt = (None, None, None)
    states[i,14] = controller_state["rw_cmd"]
    # event detection (kept as original)
    T = thrust_time(times[i]); m = mass_time(times[i]); weight = m*9.80665
    gp = np.radians(states[i,12]); gy = np.radians(states[i,13])
    tb_x = T * sin(gy); tb_y = -T * sin(gp); tb_z = T * cos(gp) * cos(gy)
    thr_world_check = body_to_world(np.array([tb_x, tb_y, tb_z]), states[i,10], states[i,6], states[i,8])
    Tz_world = thr_world_check[2]
    vz_val = states[i,5]; vz_prev = states[i-1,5]
    if flight_state == "GROUND" and Tz_world > 0.1*weight and (vz_val - vz_prev) > 0.0:
        flight_state = "BOOST"; boost_start = times[i]; events.append({"time":times[i],"event":"LAUNCH"})
    if flight_state == "BOOST" and times[i] - (boost_start or 0) >= burn_time:
        flight_state = "COAST"; events.append({"time":times[i],"event":"BURNOUT"})
    z = states[i,2]
    if z > max_alt: max_alt = z
    if flight_state == "COAST" and z < max_alt - 0.5 and apogee_time is None and max_alt > 0.5:
        flight_state = "APOGEE"; apogee_time = times[i]; events.append({"time":times[i],"event":"APOGEE"}); parachute_deploy_alt = max_alt
    if apogee_time and (not parachute_deployed) and z < parachute_deploy_alt * 0.6:
        parachute_deployed = True; params["Cd"] = 1.8; params["A"] = 0.5; flight_state = "DESCENT"; events.append({"time":times[i],"event":"PARACHUTE_DEPLOY"})
    if flight_state == "DESCENT" and z <= 0.2:
        flight_state = "LANDED"; events.append({"time":times[i],"event":"LANDED"}); states = states[:i+1,:]; times = times[:i+1]; break

# -------------------------
# POSTPROCESS: extract arrays
# -------------------------
x = states[:,0]; y = states[:,1]; alt = states[:,2]
pitch = states[:,6]; yaw = states[:,8]; roll = states[:,10]
servo_pitch = states[:,12]; servo_yaw = states[:,13]; rw_cmd = states[:,14]
thrust = np.array([thrust_time(t) for t in times])
mass = np.array([mass_time(t) for t in times])
cp_vals = np.array([cp_offset_time(t) for t in times])
ap_idx = np.argmax(alt); ap_t = times[ap_idx]; ap_alt = alt[ap_idx]

print("=== EVENT TIMELINE ===")
for e in events:
    print(f"{e['time']:.3f}s : {e['event']}")
print("=======================")
print(f"Liftoff detected? {'LAUNCH' in [ev['event'] for ev in events]}")
print(f"Max altitude: {ap_alt:.2f} m at t={ap_t:.3f}s")

# -------------------------
# Make the plots (unchanged)
# -------------------------
# plt.figure(figsize=(9,3)); plt.plot(times, thrust, label='Thrust (N)'); plt.ylabel('Thrust (N)'); plt.xlabel('Time (s)'); plt.grid(True); plt.legend(); plt.show()
# plt.figure(figsize=(9,3)); plt.plot(times, alt, label='Altitude (m)'); plt.ylabel('Altitude (m)'); plt.xlabel('Time (s)'); plt.grid(True)
# for e in events: plt.axvline(e["time"], linestyle='--', color='red'); plt.text(e["time"], max(alt)*0.9, e["event"], rotation=90)
# plt.legend(); plt.show()
# plt.figure(figsize=(9,3)); plt.plot(times, np.degrees(pitch), label='Pitch (deg)'); plt.plot(times, np.degrees(yaw), label='Yaw (deg)'); plt.ylabel('Angle (deg)'); plt.xlabel('Time (s)'); plt.grid(True); plt.legend(); plt.show()
# plt.figure(figsize=(9,3)); plt.plot(times, servo_pitch, label='Servo Pitch (deg)'); plt.plot(times, servo_yaw, label='Servo Yaw (deg)'); plt.ylabel('Servo (deg)'); plt.xlabel('Time (s)'); plt.grid(True); plt.legend(); plt.show()

# -------------------------
# CZML EXPORT (Kennedy) — no GitHub, local file only
# -------------------------
def generate_czml(times_arr, lat_arr, lon_arr, alt_arr, model_url, out_path="trajectory.czml"):
    start_time = datetime.utcnow()
    start_iso = start_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    end_iso = (start_time + timedelta(seconds=float(times_arr[-1]))).strftime("%Y-%m-%dT%H:%M:%SZ")

    czml_doc = [
        {
            "id": "document",
            "version": "1.0",
            "clock": {
                "interval": f"{start_iso}/{end_iso}",
                "currentTime": start_iso,
                "multiplier": 1,
                "range": "LOOP_STOP"
            }
        },
        {
            "id": "rocket",
            "name": "Rocket",
            "availability": f"{start_iso}/{end_iso}",
            "model": {
                "gltf": model_url,
                "scale": 0.001,  # Adjust scale as needed
                "minimumPixelSize": 64,
            },
            "position": {
                "epoch": start_iso,
                "cartographicDegrees": []
            },
            "orientation": {
                "velocityReference": "#position"
            },
            "path": {
                "material": {"solidColor": {"color": {"rgba": [255, 0, 0, 255]}}},
                "width": 3,
                "leadTime": 0,
                "trailTime": 60
            },
            "viewFrom": {
                "cartesian": [0, -500, 200]
            }
        }
    ]

    for t, lat, lon, alt in zip(times_arr, lat_arr, lon_arr, alt_arr):
        czml_doc[1]["position"]["cartographicDegrees"] += [float(t), float(lon), float(lat), float(alt)]

    with open(out_path, "w") as f:
        json.dump(czml_doc, f, indent=2)

    print(f"✅ Saved CZML to {out_path} ({len(times_arr)} samples)")

# -----------------------------
# Example: convert local XY to lat/lon around Kennedy
# -----------------------------
lat0, lon0 = 28.583434, -80.582876
deg_per_m_lat = 1.0 / 111320.0
deg_per_m_lon = 1.0 / (111320.0 * math.cos(math.radians(lat0)))

lats = lat0 + (y * deg_per_m_lat)
lons = lon0 + (x * deg_per_m_lon)
alts_m = alt

# Local rocket GLTF file path (relative to viewer.html)
rocket_model = "models/Rocket.glb"

# Save CZML locally (no git, no remote)
generate_czml(times, lats, lons, alts_m, rocket_model, out_path="trajectory.czml")

# -----------------------------
# Create viewer.html automatically (only once if not exists)
# -----------------------------
# Load Kennedy trajectory CZML from file
import http.server
import socketserver
import webbrowser
import os

# Path to your Kennedy trajectory CZML file
CZML_FILE = "trajectory.czml"  # Make sure this is in the same folder

# Read the Kennedy trajectory CZML
with open(CZML_FILE, "r", encoding="utf-8") as f:
    trajectory_czml = f.read()

# HTML content with rocket model + Kennedy trajectory + guaranteed globe imagery
html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Kennedy Rocket Launch</title>
    <script src="https://cesium.com/downloads/cesiumjs/releases/1.119/Build/Cesium/Cesium.js"></script>
    <link href="https://cesium.com/downloads/cesiumjs/releases/1.119/Build/Cesium/Widgets/widgets.css" rel="stylesheet">
    <style>
        html, body, #cesiumContainer {{ width: 100%; height: 100%; margin: 0; padding: 0; }}
    </style>
</head>
<body>
<div id="cesiumContainer"></div>
<script>
    Cesium.Ion.defaultAccessToken = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJqdGkiOiI0YTg5MDIxZS01N2I0LTRmYjItOGRkZC00MTI3ZWM3MzJlNzYiLCJpZCI6MzMxMzM3LCJpYXQiOjE3NTUxODA5Nzh9.5mG6UMABdHIBXtrxRR5554zR2v-PSAquNtPt4zetNvk";

    var viewer = new Cesium.Viewer('cesiumContainer', {{
        imageryProvider: Cesium.createWorldImageryAsync({{
            style: Cesium.IonWorldImageryStyle.AERIAL
        }}),
        shouldAnimate: false
    }});

    var czmlData = {trajectory_czml};

    // Load the Kennedy trajectory
    viewer.dataSources.add(Cesium.CzmlDataSource.load(czmlData)).then(function(ds) {{
        viewer.trackedEntity = ds.entities.values[0]; // Track first entity from file
    }});
</script>
</body>
</html>
"""

# Save HTML file
with open("index.html", "w", encoding="utf-8") as f:
    f.write(html_content)

# Start server and open browser
PORT = 8000
Handler = http.server.SimpleHTTPRequestHandler

with socketserver.TCPServer(("", PORT), Handler) as httpd:
    print(f"Serving at http://localhost:{PORT}")
    webbrowser.open(f"http://localhost:{PORT}/index.html")
    httpd.serve_forever()
