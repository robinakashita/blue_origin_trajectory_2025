# Blue Origin Suborbital Trajectory Simulation (Physics 360 Final Project)

This repository implements a numerical simulation of a Blue Origin–style suborbital ascent, coast, and descent trajectory, calibrated using publicly available parameters from the New Shepard flight on April 14, 2025. The simulation solves the coupled equations of motion for a vertically launched rocket with variable mass, atmospheric drag, and altitude-dependent gravity, and reproduces key characteristics of an ~11-minute suborbital flight to the Karman line.

The model is implemented in a single, self-contained Jupyter Notebook:

- `Final_Report_Robin_N.ipynb`

---

## 1. Scientific Motivation

Suborbital launch vehicles provide an accessible test bed for studying rocket dynamics, aerodynamic drag, and ballistic flight under realistic atmospheric conditions. Although Blue Origin’s exact parameters are proprietary, approximate public values (mass, thrust, burnout altitude, etc.) are sufficient to reproduce a physically reasonable and safe trajectory.

This project demonstrates how undergraduate-level physics and numerical methods can be combined to:

- Construct a time-varying force balance for a multistage rocket,
- Quantify maximum altitude, flight duration, and time spent in space, and
- Analyze g-loads experienced during ascent and re-entry.

---

## 2. Physical Model

The rocket is modeled as a point mass moving in the vertical direction \( y(t) \). The governing equation of motion is

$$
m(t)\frac{dv}{dt} = T(t) - D(v,h) - m(t)g(h)
$$

where $m_{\text{prop}}(t)$ decreases linearly during the burn phase and is constant thereafter (booster separation / burnout).

- $( m(t) )$ is the time-dependent mass (structure + propellant),
- $( T(t) )$ is the thrust (nonzero during powered ascent, zero during coast and descent),
- $( D(v, h) )$ is the aerodynamic drag force, and
- $( g(h) )$ is the gravitational acceleration as a function of altitude \( h \).

### 2.1 Mass Model

The total mass is decomposed as

$$
m(t) = m_{\text{dry}} + m_{\text{prop}}(t),
$$

where $m_{\text{prop}}(t)$ decreases linearly during the burn phase and is constant thereafter (booster separation / burnout).

### 2.2 Thrust

Thrust is approximated as a piecewise function

$$
T(t) =
\begin{cases}
T_0, & 0 \le t \le t_{\text{burn}} \\
0,   & t > t_{\text{burn}},
\end{cases}
$$

with parameters chosen to yield realistic burnout altitude and velocity.

### 2.3 Drag

Aerodynamic drag is modeled using the standard quadratic form

$$
D(v, h) = \frac{1}{2}\rho(h)v^2C_dA,
$$

where

- $( C_d )$ is the drag coefficient,
- $( A )$ is the reference cross-sectional area, and
- $( \rho(h) )$ is the atmospheric density as a function of altitude.

### 2.4 Atmospheric Density and Gravity

Atmospheric density decreases exponentially with height,

$$
\rho(h) = \rho_0 \exp\left(-\frac{h}{H}\right),
$$

where $( \rho(h) )$ is sea-level density and $( H )$ is the scale height.

Gravitational acceleration is modeled using an inverse-square dependence:

$$
g(h) = g_0 \left(\frac{R_E}{R_E + h}\right)^2,
$$

with $R_E$ the Earth’s mean radius and $g_0 \approx 9.81\mathrm{m\/s^{-2}}$.

---

## 3. Numerical Methods

The equations of motion are integrated using `scipy.integrate.solve_ivp` over three main phases:

1. Powered ascent (engine on, mass decreasing),
2. Coast and ballistic ascent (engine off, drag + gravity only),
3. Descent and re-entry (including increased drag and parachute deployment if modeled).

The notebook computes and stores:

- Time array $( t )$,
- Altitude $( y(t) )$,
- Horizontal displacement $( x(t) )$ (nominally near zero for a vertical trajectory),
- Velocity components $( v_x(t), v_y(t) )$,
- Acceleration magnitude $( |a(t)| )$.

From these, the code derives:

- Maximum altitude),
- Maximum speed,
- Booster separation time and altitude,
- Time spent above the Karman line ($h \ge 100\\mathrm{km}$),
- Total flight duration.

---

## 4. Key Results

Using the current parameter set, the simulation yields:

- Maximum altitude : ~137 km  
- Booster separation: ~141 s after launch at ~51.9 km altitude  
- Maximum ascent speed: ~1,360 m/s  
- Total flight duration: ~647 s (~10.8 minutes)  
- Time above 100 km (Karman line): ~177 s (~3 minutes)  

These values are broadly consistent with public descriptions of New Shepard flights and illustrate the physical regimes encountered in suborbital trajectories.

---

## 5. Figures

- **Position vs Time**
  
![Position vs Time](figures/position_vs_time.png)

- **Velocity vs Time**

![Velociy vs Time](figures/velocity_vs_time.png)


- **Acceleration vs Time**

![Acceleration vs Time](figures/acceleration_vs_time.png)
