"""
@file   mudkip.py

Based on input a given field configuration and hourly weather data, optimize the
input of water by an irrigation system such that the available groundwater 
volume remains greater than a given level throughout the day.

@author Jake Lee/nubby  (jlee211@ucsc.edu)
@date   13 Dec 2024
"""
import csv
import sys

import cvxpy as cp
import numpy as np

import matplotlib.pyplot as plt


# Module-level variables.
_DEFAULT_FARM_DATA_PATH = "./data/sample_farm_data.csv"
_STEEP_FARM_DATA_PATH = "./data/sample_steep_data.csv"
_DEFAULT_WEATHER_DATA_PATH = "./data/sample_weather_data.csv"
_DEFAULT_SCHEDULE_PATH = "./data/sample_flat_schedule.csv"
_STEEP_SCHEDULE_PATH = "./data/sample_steep_schedule.csv"
_DEFAULT_DURATION = 24   # Hours

# Datatypes.
class Farm(object):
    """
    Each Farm can be represented by an n-sized array of constituent Field nodes,
    F, and a corresponding diagonally-symmetrical matrix of size n x n, A. Each
    element, A[i][j], describes the flow (change in water volume) that occurs
    between node i and j. These values are derived from differences in
    hydraulic potential between nodes.

    The interactions between each node can be modeled by the Criteria3D model
    of water dynamics, namely:

    dV_i/dt = sum(Q_ij + q_i) for all i != j;
    Q_ij = -K_ij * S_ij * (H_i - H_j) / L_ij.

    In the above, "q_i" is the input water at a given node, which is what we
    aim to minimize (subject to some minimum soil water volume being preserved
    at each node).

    The rest of the parameters:
    + V_i:  Amount of stored water in volume surrounding node i.
    + Q_ij: Flux between the ith and jth node.
    + K_ij: Harmonic mean of the hydraulic conductivity at the ith and jth
            nodes.
    + S_ij: The interfacial area between nodes i and j; the area of the plane
            between the surface and the groundwater basin that can be drawn
            between nodes i and j.
    + L_ij: The distance between nodes i and j.

    @param  data    (list[dict])  Farm data [that we care about for now] in dict
                                  format:
        [
            {
                "Name": str,
                "PoreSizeDistFactor": float,
                "RetentionInflection": float,
                "SaturatedHydraulicConductivity": float,
                "SW": float,    <-- "SurfaceWaterVolume"
                "X": float,
                "Y": float,
                "Z": float      <-- "Elevation"
            }, ...
        ]
    """
    def __init__(self, data: list[dict]):
        self.F = [] 
        self.Q = None
        self.size = len(data)
        self._build_farm(data)

    def __repr__(self):
        """
        """
        return (
            str([str(field.Z) for field in self.F])
        )

    def _build_farm(self, data: list[dict]):
        """_build_farm(size)
        Initialize a new array of Fields of a desired size and build their
        adjacency matrix.
        """
        print("Building farm...")
        # Generate and initialize each Field node.
        [self.F.append(Field(field_data)) for field_data in data]
        # Create an associative, diagonally-symmetric matrix.
        n = len(data)
        self.Q = np.zeros((n,n), dtype="float")
        self.update()
        print("    DONE.")

    def _calculate_Q_ij(self, i: int, j: int):
        """_calculate_Q_ij(i,j)
        Calculate the change in water volume at node i due to node j.
        """
        Fi = self.F[i]
        Fj = self.F[j]
        
        # Distances. <-- For very large farms, we can save these here.
        Fi_coords = np.array([Fi.X, Fi.Y, Fi.Z])
        Fj_coords = np.array([Fj.X, Fj.Y, Fj.Z])
        L_ij = np.linalg.norm(Fi_coords - Fj_coords)

        # Harmonic mean of Ki and Kj.
        dK = 2 / ((1 / Fi.K) + (1 / Fj.K))

        # Water interfacial area between nodes.
        Fi_coords = np.array([Fi.X, Fi.Y])
        Fj_coords = np.array([Fj.X, Fj.Y])
        d_horiz = np.linalg.norm(Fi_coords - Fj_coords)
        S_ij = min([Fi.S, Fj.S]) * d_horiz
        if (Fi.S != Fj.S):
            S_ij += (1/2) * abs(Fi.S - Fj.S) * d_horiz
        
        # Return Q_ij.
        Q_ij = -dK * S_ij * (Fi.H - Fj.H) / L_ij
        return Q_ij

    def update(self):
        """update()
        Update the representative matrix with the status of every constituent
        Field node.
        """
        n = len(self.F)
        for i in range(n):
            for j in range(n):
                if (i == j):
                    self.Q[i, j] = 0
                else:
                    self.Q[i, j] = self.Q[i, j] + self._calculate_Q_ij(i, j)


class Field(object):
    """Field()
    Container for information and dynamics relating to each Field node.
    NOTE: These aim to include, if not add to, all of the properties of Field
          nodes used in APSIM.

    A Field imports the following data:

        {
            "Name": str,
            "PoreSizeDistFactor": float,
            "RetentionInflection": float,
            "SaturatedHydraulicConductivity": float,
            "SW": float,    <-- "SurfaceWaterVolume"
            "X": float,
            "Y": float,
            "Z": float      <-- "Elevation"
        }
    """
    def __init__(self, data: dict):
        try:
            self.Name = data["Name"]
            self.SurfaceWaterVolume = float(data["SW"])
            # Get position data.
            self.X = float(data["X"])
            self.Y = float(data["Y"])
            self.Z = float(data["Z"])
            # Set up field params.
            self.K = 1
            self.H = 0
            self.S = 1    # Surface water depth at node.
            self._build_field(data)
        except KeyError as e:
            raise(e)

    def __repr__(self):
        return "{}: Elevation = {}m; Water Volume = {}".format(
            self.Name,
            self.Z,
            self.SurfaceWaterVolume
        )

    def _build_field(self, data: dict):
        """_build_field(data)
        Generate Field parameters, including:
        + H =   Total hydraulic head at node, given by the sum of the
                gravitational and hydraulic matric potentials.
        + K =   Hydrolic conductivity at node.
        """
        # Set up static variables for calculating Field params.
        self.alpha = float(data["RetentionInflection"])
        self.h_g = self.Z       # Gravitational head is simply the elevation.
        self.g = 9.8            # m/s^2
        self.K_s = float(data["SaturatedHydraulicConductivity"])
        self.l = 0.5            # From CRITERIA3D.
        self.n = float(data["PoreSizeDistFactor"]) # Unitless.
        self.m = 1 - 1 / self.n # From CRITERIA3D.

        self._update_H()
        self._update_K()
        self._update_S()

    def _update_H(
            self,
            soil_h2o_density: float = 0.0015,
            soil_matric_potential: float = 0.5
        ):
        """_update_H(soil_h2o_density, soil_matric_potential)
        Update the current node's hydraulic head.
        NOTE: Currently assumes that there is no standing water at node.
        
        @param soil_h2o_density         (float)
        @param soil_matric_potential    (float)
        """
        ro = soil_h2o_density           # kg/m^3
        # TODO(nubby)
        p_ms = soil_matric_potential    # m
        # Pressure head initializes at 0, but is determined by soil water
        # so we should recalculate H each time step.
        h_p = p_ms / (self.g * ro)
        self.H = h_p + self.h_g

    def _update_K(self):
        """
        Update the current node's hydraulic conductivity based on saturated
        conductivity, water retention parameters, and water potential.
        """
        self.K = self.K_s * (
            1 - (self.alpha * self.H)**(self.m * self.n) * (
                1 + (self.alpha * self.H)**self.n
            )**(-self.m)
        ) / (1 + (self.alpha * self.H)**(self.n))**(self.m * self.l)

    def _update_S(self):
        """
        Update the current node's water flow height.
        NOTE: For now, this value is kept low and constant. The literature
              defines this parameter as a coupling of surface water flow depth
              and the gradual flow of water through wet soil beneath the
              surface. On a hot day on a farm, there will be no water on the
              surface, but the soil should stay moist; hence, this needs to have
              at least a small value. Also, without a small value, all water
              flow between nodes immediately becomes 0. <-- TODO
        """
        # Pretty sure S is the surface area of the Field.
        self.S = 1

# Initialize H and K.  self._update_H() self._update_K()
    def irrigate(self, duration: float, rate: float):
        """irrigate(duration, rate)
        Models BOTH rainfall and user input irrigation.
        Updates self.S, which indicates surface water depth in m.
        """
        pass


# Helpers.
def _load_csv(fpath: str) -> list[dict]:
    data = []
    print(f"Reading from {fpath}...")
    with open(fpath, "r+") as csvs:
        reader = csv.DictReader(csvs)
        for row in reader:
            data.append(row)
    print(f"    DONE")
    if (not data):
        raise OSError(
            f"WARNING!! {fpath} is an empty file!"
        )
    return data

def _write_csv(fpath: str, data: list):
    print(f"Writing to {fpath}...")
    with open(fpath, "w+") as csvs:
        writer = csv.writer(csvs, delimiter=",")
        for row in data:
            writer.writerow(row)
    print(f"    DONE")


# Main functionality.
def FAO(hour: int, T: float) -> float:
    """FAO
    Convert a temperature (C) into an evapotranspiration rate using the FAO
    Penman-Monteith equation.

    @param hour (int)   Integer hour closest to the time of measurement.
    @param T    (float) Temperature in C at that time.
    """
    e_o = 0.6108 * np.exp(17.27 * T / (T + 237.3))
    delta = 4098 * e_o / (T + 237.3)**2
    R = (-(T - 12)**2 + 40) / 24
    psych = 0.665 * 10**(-3)
    # Assuming no wind speed for this simulation.
    u_w = 0
    # Assuming there is no saturation pressure deficit.
    e_s = 0
    e_a = 0
    # Assuming that there is no soil heat flux.
    G = 0

    # Calculate ET rate.
    et = (0.408 * delta * (R - G) + psych * (
        900 / (T + 273)) * u_w * (e_s - e_a)
    ) / (delta + psych * (1 + 0.34 * u_w))
    return et

def _subplot_config(f, i: int, n: int, u, V_g, V_s):
    """_subplot_config(n, i, u, V_g, V_s)
    """
    index = i - 1
    ax1 = f.add_subplot(n, 1, i)
    ax1.plot(
        V_g[index, :].value + V_s[index, :].value, c="c"
    )
    ax1.set_ylabel(f"$V(t)_{i}$", fontsize=16, c="c")
    ax1.set_xticks([])
    ax2 = ax1.twinx()
    ax2.plot(u[index, :].value, c="r")
    ax2.set_ylabel(f"$u(t)_{i}$", fontsize=16, c="r")
    return ax1
    

def plot_results(n, u, V_g, V_s):
    """plot_results(n, u, V)
    """
    f = plt.figure(figsize=(10,8))

    ax1 = _subplot_config(f, 1, n, u, V_g, V_s)
    ax1.set_title("V (m^3) and u (m^3) vs Hour of Day", fontsize=20)
    
    [_subplot_config(f, i, n, u, V_g, V_s) for i in range(2,n)]

    ax1 = _subplot_config(f, n, n, u, V_g, V_s)
    ax1.set_xticks(range(0,len(u.value[0])+1))

    plt.show()

def _format_schedule(data) -> dict[list[str]]:
    time = lambda t: str(t).rjust(2,"0")+":00"
    schedule = []
    header = [
        "Time",
        "Field1",
        "Field2",
        "Field3",
        "Field4",
        "Field5",
        "Field6",
        "Field7",
        "Field8",
        "Field9"
    ]
    schedule.append(header)
    for t in range(len(data[0])):
        schedule.append([
            time(t),
            round(data[0][t], 2),
            round(data[1][t], 2),
            round(data[2][t], 2),
            round(data[3][t], 2),
            round(data[4][t], 2),
            round(data[5][t], 2),
            round(data[6][t], 2),
            round(data[7][t], 2),
            round(data[8][t], 2)
        ]) 

    return schedule


def run_mudkip(farm: Farm, weather: list[dict]):
    """run_mudkip(farm)
    When is the best time to water your crops?...on a hill.

    @param  farm        (Farm)  The structure/dynamics of the simulation.
    @return schedule    (dict)  
    """
    n = len(farm.F)
    Q = farm.Q
    T = _DEFAULT_DURATION
    tdelta = T / 24

    et = []

    # Threshold volume needed at each Field; doubles as initial value.
    alpha = np.full(n, 1, dtype="float")
    # u = input water at each Field node.
    u = cp.Variable((n, T))
    u_0 = np.full(n, 0, dtype="float") # Start without any irrigation for now.
    u_max = np.full(n, 25, dtype="float")
    et_0 = np.full(n, FAO(0, float(weather[0]["Temperature"])), dtype="float")
    q_0 = u_0 + et_0

    # These parameters simulate the relative moisture retention factors of air
    # and soil.
    sigma = 0.3     # Surface water evaporates quickly.
    epsilon = 0.6   # Water stays longer in soil.
 
    V_s0 = np.full(n, 1)
    V_g0 = np.full(n, 2)

    V = cp.Variable((n, T + 1))
    V_s = cp.Variable((n, T + 1))
    V_g = cp.Variable((n, T + 1))
    
    cost = 0
    constr = []

    # Calculate the sum(Q_i) array to describe water dynamics for each node.
    sQ = []
    for i in range(n):
        Q_i = np.sum(Q[i])
        sQ.append(Q_i)

    # Add a water storage parameter to this to account for soil absorption.
    # TODO(nubby): Bind to characteristics of each Field node type.
    S = np.full(n, 0)

    for t in range(T):
        # Calculate the rate of evapotranspiration at t.
        et =  FAO(t, float(weather[t]["Temperature"]))
        constr += [
            V[:, t+1] == V_s[:, t+1] + V_g[:, t+1],
            V_g[:, t+1] == V_g[:, t] + epsilon * V_s[:, t] + (1-epsilon) * et,
            V_s[:, t+1] == (1-epsilon) * V_s[:, t] + (1-sigma) * et + (
                u[:, t] + Q @ V_s[:, t]
            )
        ]
        constr += [V_g[:, t] >= alpha]
        # Add the sum of the water inputs to the cost at each step.
        cost += cp.sum(u[:, t])
        # TODO(nubby): The following constraints make the primal problem
        #              unsolveable; try to reframe this using the Lagrangian
        #              dual.
        #constr += [u[:, t] >= 0]
        constr += [u[:, t] >= 0, u[:, t] <= u_max]
        #constr += [V_s[:, t] >= 0]
    constr += [V[:, 0] == V_s0 + V_g0]
    constr += [V_g[:, 0] == V_g0, V_s[:, 0] == V_s0]
    constr += [V_g[:, T] == V_g0]
    constr += [u[:, 0] == u_0]

    problem = cp.Problem(cp.Minimize(cost), constr)
    problem.solve()

    print("Problem status: ", problem.status)
    print("Total water used: ", problem.value)
    #print("Water volumes: ", V_g.value + V_s.value)
    print("Watering schedule: ", str(u.value))

    # Reporting to STDOUT.
    av = 0
    for i, s in enumerate(u.value):
        print(i,". ",sum(s));
        av += sum(s)
    av /= len(u.value)
    print("Average (m^3): ",av)

    plot_results(n, u, V_g, V_s)

    results = _format_schedule(u.value)

    return results 

def clean_up(e: Exception):
    print(f"{e}; good bye.")
    sys.exit(0)

if __name__ == "__main__":
    # Load location-identifiable data corresponding to a farm being studied.
    try:
        farm_data = _load_csv(_DEFAULT_FARM_DATA_PATH)
        #farm_data = _load_csv(_STEEP_FARM_DATA_PATH)
    except OSError as e:
        clean_up(e)

    # Generate a model of the farm being explored.
    try:
        farm = Farm(farm_data)
    except ValueError or KeyError as e:
        clean_up(e)

    # Load a weather "schedule", i.e. the number of time steps and amount of
    # rain to simulate on the farm.
    try:
        weather_data = _load_csv(_DEFAULT_WEATHER_DATA_PATH)
    except OSError as e:
        clean_up(e)

    # Generate an optimal irrigation schedule for each Field in the simulation.
    watering_schedule = run_mudkip(farm, weather_data)

    # Output a watering schedule for each node.
    # NOTE: Watering "schedules" just describe the total volume of water that
    #       (m^3) should be output by each sprinkler during the course of a
    #       given hour to achieve optimal water efficiency for each node.
    _write_csv(_DEFAULT_SCHEDULE_PATH, watering_schedule)
    #_write_csv(_STEEP_SCHEDULE_PATH, watering_schedule)
