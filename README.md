# Modeling the Ascent and Descent of Blue Origin’s New Shepard Capsule Using Numerical Methods in 2-D (Physics 360 Final Project)

This project aims to analyze the trajectory of the Blue Origin's New Shepard capsule's 11-minute journey from the ground to the Karman Line and back to the ground. This spacecraft will experience a point of maximum stress, main engine cut-off, capsule separating from the booster, capsule passing the Karman Line, parachutes deploying, and finally, the capsule touching down on the landing pad. This system does not have an analytic solution for all conditions, which makes it a great idea for a computational physics project. I will use numerical integration methods (Runge-Kutta) to solve the equations of motion and (Scipy/scipy.integrate.solve_ivp) to calculate the time step for each step. The model will simulate position, velocity, and acceleration as a function of time, and it will estimate how long the space crew experienced weightlessness after the capsule passed the Karman Line.

The model is implemented in a single, self-contained Jupyter Notebook:

- `Physics360_Final_Report_Robin_N.ipynb`

---

## 1. Scientific Motivation

As space travel becomes more accessible to the general public, it becomes crucial to study the effects of Newton's Second Laws of Motion and enhance the safety of spacecraft trajectories and space missions. By modeling the trajectory of Blue Origin's New Shepard Capsule from the beginning to the end, we are able to visualize where the maximum terminal velocity of the capsule occurs. By splitting up our 2nd-Order ODE into two First-Order ODEs, we are able to numerically calculate the maximum terminal velocity for the Blue Origin capsule and the free-fall time for the space crew. The capsule will travel from $t = 0 s$ to $t = 660s$. At $t = 141s$, the engine of the capsule will cut-off, and the booster will separate from the capsule. The position and velocity will start at $x$, $y$, $v_x$, $v_y$ = $0$, respectively. The mass of the capsule will start at $m = 75,000 kg$.

Sources:
- [Blue Origin For the Benefit of Earth](https://www.blueorigin.com/new-shepard)
- [SciPy](https://scipy.org/) and Matplotlib documentation  – I am planning to use Python with these libraries for numerical integration and plotting.
- [Runge-Kutta method](https://math.libretexts.org/Courses/Monroe_Community_College/MTH_225_Differential_Equations/03%3A_Numerical_Methods/3.03%3A_The_Runge-Kutta_Method) will be implemented to solve the differential equation.
- New Shepard Flight Test Results from Blue Origin De-Orbit Descent and Landing Tipping Point for spacecraft (https://doi.org/10.2514/6.2022-1829)

This project demonstrates how undergraduate-level physics and numerical methods can be combined to:

- Construct a time-varying force balance for a multistage rocket,
- Quantify maximum altitude, flight duration, and time spent in space, and
- Analyze g-loads experienced during ascent and re-entry.

---

## 2. Physical Model

The rocket is modeled as a point mass moving in the horizontal and vertical direction \( x(t) and y(t) \). The governing equation of motion is

Newton's Second Law:

$$
m_{\rm cap}\\frac{d^2 \vec{r}}{dt^2}
\=\
\underbrace{F(t)\\hat{\jmath}}_{\displaystyle\mathbf F_{\rm thrust}}
\+\
\underbrace{\Bigl(-\tfrac12\\rho(h)C_DA\|\mathbf v|\\mathbf v\Bigr)}_{\mathbf F_{\rm drag}}
\+\
\underbrace{\Bigl(+\tfrac12\\rho(h)C_LA\|\mathbf v|^2\\hat{\jmath}\Bigr)}_{\mathbf F_{\rm lift}}
\+\
\underbrace{\Bigl(0, - m_{\rm cap}\ g(h)\Bigr)}_{\mathbf{F}_{\rm grav}}.
$$

$$
m_{\text{cap}} \frac{d^2 \vec{r}}{dt^2} = \vec{F}_{\text{thrust}} + G \cdot \frac{m_{\text{cap}}\,m_E}{(R_E + h)^2} - \frac{1}{2} \rho(h) C A v^2 + \vec{F}_{\text{lift}}
$$

$$ \vec{r}= \sqrt{x^2 + y^2}\ $$

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

## 3. Define All Parameters

### Define All Parameters:
#### 1. State Vector

$$
\vec{Y}(t) = 
\begin{bmatrix}
x(t) \\
y(t) \\
v_x(t) \\
v_y(t) \\
m
\end{bmatrix}
$$

#### 2. 2D Equations of Motion with Drag and Gravity
##### Velocity in X-Direction:
$$ 
\dot{x} = v_x(t) 
$$
##### Velocity in Y-Direction:
$$ 
\dot{y} = v_y(t) 
$$
##### Mass Flow Rate
$$ 
\dot{m} = -\frac{F(t)}{I_{Isp}g_0} 
$$ 

#### 3. Split Up the 2nd-Order ODE into two  First Order ODE's
##### Acceleration in X-Direction:
$$ 
\dot{v_x} = -\frac{1}{2m_{\rm cap}}\\rho(y)C_DAvv_x 
$$ 
##### Acceleration in Y-Direction:
$$ \dot{v_y} = \frac{F(t)}{m_{\rm cap}}
   \-g(y)
   \-\\frac{1}{2m_{\rm cap}}\\rho(y)C_DAvv_y
   \+\\frac{1}{2m_{\rm cap}}\\rho(y)C_LAv^2 
$$
   
---

## 4. Numerical Methods

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

- Maximum altitude,
- Maximum speed,
- Booster separation time and altitude,
- Time spent above the Karman line ($h \ge 100\\mathrm{km}$),
- Total flight duration.

---

## 5. Key Results

Using the current parameter set, the simulation yields:

- Maximum altitude : ~137 km  
- Booster separation: ~141 s after launch at ~51.9 km altitude  
- Maximum ascent speed: ~1,360 m/s  
- Total flight duration: ~647 s (~10.8 minutes)  
- Time above 100 km (Karman line): ~177 s (~3 minutes)  

These values are broadly consistent with public descriptions of New Shepard flights and illustrate the physical regimes encountered in suborbital trajectories.

---

## 6. Figures

- **Position vs Time**
  
![Position vs Time](figures/position_vs_time.png)

- **Velocity vs Time**

![Velociy vs Time](figures/velocity_vs_time.png)


- **Acceleration vs Time**

![Acceleration vs Time](figures/acceleration_vs_time.png)
