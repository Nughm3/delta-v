from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from scipy.optimize import minimize
from sgp4.api import Satrec
from functools import partial

import math
import numpy as np

INPUT = "data/cosmos-1408.tle"
OUTPUT = "delta-v/cosmos-1408.csv"

EPOCHS = 50 # days
TRANSFER_MAX = 5 # days

# SI units
J2 = 1.08262668e-3
G = 6.6743015e-11
M = 5.972e24
R = 6378137


@dataclass(frozen=True)
class Debris:
    sma: float  # m
    incl: float  # radians
    raan: float  # radians

    @staticmethod
    def from_tle(tle):
        if len(tle) == 3:
            s, t = tle[1:]
        else:
            s, t = tle

        sat = Satrec.twoline2rv(s, t)
        return Debris(sat.a * sat.radiusearthkm * 1000, sat.inclo, sat.nodeo)


def raan_constraint(i, j, k, m):
    SECONDS_IN_DAY = 60 * 60 * 24
    k *= SECONDS_IN_DAY
    m *= SECONDS_IN_DAY

    def constraint_function(x):
        sma, incl = x[0], x[1]
        raan_chaser = (
            i.raan
            + nodal_precession(i.sma, i.incl) * k
            + nodal_precession(sma, incl) * m
        ) % math.tau
        raan_j = (j.raan + nodal_precession(j.sma, j.incl) * (k + m)) % math.tau

        return raan_chaser - raan_j

    return constraint_function


def nodal_precession(sma, incl):
    global J2, G, M, R
    return -3 / 2 * J2 * math.sqrt(G * M) * R**2 * sma ** (-7 / 2) * math.cos(incl)


def delta_v(x, i, j):
    global G, M

    sma, incl = x[0], x[1]

    di_i = abs(i.incl - incl)
    di_j = abs(j.incl - incl)

    s_i = (
        0
        if di_i == 0
        else (1 / di_i)
        * math.atan(math.sin(di_i) / (math.sqrt(sma**3 / i.sma**3) + math.cos(di_i)))
    )
    s_j = (
        0
        if di_j == 0
        else (1 / di_j)
        * math.atan(math.sin(di_j) / (math.sqrt(sma**3 / j.sma**3) + math.cos(di_j)))
    )

    v_i = math.sqrt(G * M / i.sma)
    v_j = math.sqrt(G * M / j.sma)
    v_d = math.sqrt(G * M / sma)

    sma_id = (i.sma + sma) / 2
    sma_dj = (sma + j.sma) / 2

    vh1p = vis_viva(i.sma, sma_id)
    vh1a = vis_viva(sma, sma_id)
    vh2p = vis_viva(sma, sma_dj)
    vh2a = vis_viva(j.sma, sma_dj)

    dV_h1p = cosine_rule(v_i, vh1p, di_i * s_i)
    dV_h1a = cosine_rule(v_d, vh1a, di_i * (1 - s_i))
    dV_h2p = cosine_rule(v_j, vh2p, di_j * s_j)
    dV_h2a = cosine_rule(v_d, vh2a, di_j * (1 - s_j))

    return dV_h1p + dV_h1a + dV_h2p + dV_h2a


def vis_viva(r, a):
    global G, M
    return math.sqrt(G * M * (2 / r - 1 / a))


def cosine_rule(a, b, C):
    return math.sqrt(a**2 + b**2 - 2 * a * b * math.cos(C))


def optimal_transfer(i, j, k, m):
    x0 = np.array([(i.sma + j.sma) / 2, (i.incl + j.incl) / 2])
    constraint = {"type": "eq", "fun": raan_constraint(i, j, k, m)}
    bounds = [(6000000, None), (0, math.tau)]

    opt = minimize(
        delta_v,
        x0,
        args=(i, j),
        method="SLSQP",
        constraints=[constraint],
        bounds=bounds,
    )

    if not opt.success:
        print("error:", res.message)
        return None
    sma, incl, cost = float(opt.x[0]), float(opt.x[1]), opt.fun
    return sma, incl, cost


def calculate_optimal_transfer(debris, delta_v, i, j, k, m):
    result = optimal_transfer(debris[i], debris[j], k, m)
    if result is not None:
        sma, incl, cost = result
        delta_v[(i + 1, j + 1, k, m)] = (sma, incl, cost)
        print(f"i={i} j={j} k={k} m={m} : sma={sma}, incl={incl}, cost={cost}", end="\r")


def main():
    global INPUT, OUTPUT

    debris = []
    with open(INPUT) as f:
        data = f.read().splitlines()
        for i in range(0, len(data), 3):
            debris.append(Debris.from_tle(data[i : i + 3]))
    print(f"input read from {INPUT}")

    delta_v = {}
    n = len(debris)
    with ThreadPoolExecutor() as pool:
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                for k in range(1, EPOCHS + 1):
                    for m in range(1, TRANSFER_MAX + 1):
                        if k + m > EPOCHS:
                            break
                        pool.submit(
                            partial(calculate_optimal_transfer, debris, delta_v), i, j, k, m
                        )

    with open(OUTPUT, "w") as f:
        for k, v in delta_v.items():
            line = ",".join(list(k) + list(v))
            f.write(f"{line}\n")
    print(f"output written to {OUTPUT}")


if __name__ == "__main__":
    main()
