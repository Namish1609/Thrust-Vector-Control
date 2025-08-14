# 🚀 Model Rocket Project — Real & Simulated Flight

This repository contains:

1. **Teensy firmware** for a real model rocket (Arduino C++)  
2. **Python + CesiumJS** 3D simulation for visualizing rocket trajectories  

---

## 📂 Repository Structure

/Teensy → Teensy firmware (Arduino code)
/sim → Python scripts + CZML generator + CesiumJS viewer
/models → 3D rocket models (GLTF/GLB)
/data → Example telemetry + CZML files

yaml
Copy
Edit

---

## 🖥 Simulation Overview

The simulation:

- Generates a **CZML** file for the rocket’s trajectory  
- Runs a **local web server** to view the flight in CesiumJS  
- Displays a **3D rocket model** moving over realistic globe imagery  
- Allows pausing/resuming with Cesium’s built-in controls  

---

## 📦 Dependencies

### 🐍 Python Simulation (install via `pip`)
- `numpy`
- `matplotlib`
- `scipy`
- `requests`

Example install command:
```bash
pip install numpy matplotlib scipy requests
🔌 Teensy Firmware (install via Arduino Library Manager)
Wire (built-in)

SPI (built-in)

Servo (built-in)

SD (built-in)

Adafruit BMP3XX

MPU9250

🚀 Running the Simulation
Install the Python dependencies:

bash
Copy
Edit
pip install numpy matplotlib scipy requests
Run the main simulation script:

bash
Copy
Edit
python sim/main.py
Your browser will open showing the Cesium 3D view of the rocket flight.

📜 License
MIT — free to use and modify.

🤝 Contributions
Pull requests are welcome.
Flight data and experiments can be shared in GitHub Discussions.
