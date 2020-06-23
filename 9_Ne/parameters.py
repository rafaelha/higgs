import numpy as np

params_ = [
    {
        "Ne": np.arange(100,1400,50),
        "tmin": [-4*(2*np.pi)],
        "tmax": [25*(2*np.pi)],
        "Nt": 800,
        "T": [0.02],
        "wd":  [10],
        "s": np.array([1]),
        "m": [ np.array([1.0]) ],
        "ef": [ np.array([500]) ],
        "g": [ np.array([10]) ],
        "U": [ np.array([[0.21]]) ],
        "A0": 0.01,
        "tau": [0.2*np.pi],
        "w":  [0],
        "A0_pr": [0.01],
        "tau_pr": [2*np.pi],
        "w_pr": [1.026],
        "t_delay": [ 10*np.pi ],
    },
    {
        "Ne": np.arange(100,1400,50),
        "tmin": [-4*(2*np.pi)],
        "tmax": [25*(2*np.pi)],
        "Nt": 800,
        "T": [0.06],
        "wd":  [10],
        "s": np.array([1]),
        "m": [ np.array([1.0]) ],
        "ef": [ np.array([500]) ],
        "g": [ np.array([10]) ],
        "U": [ np.array([[0.27]]) ],
        "A0": 0.01,
        "tau": [0.3*np.pi],
        "w":  [0],
        "A0_pr": [0.01],
        "tau_pr": [0.6*np.pi],
        "w_pr": [1.026],
        "t_delay": [ 10*np.pi ],
    }
]

params = []
for p in params_:
    for A0_pr in p["A0_pr"]:
        for Ne in p["Ne"]:
            for tmin in p["tmin"]:
                for tmax in p["tmax"]:
                    for T in p["T"]:
                        for wd in p["wd"]:
                            for m in p["m"]:
                                for ef in p["ef"]:
                                    for U in p["U"]:
                                        for tau in p["tau"]:
                                            for w in p["w"]:
                                                for tau_pr in p["tau_pr"]:
                                                    for w_pr in p["w_pr"]:
                                                        for g in p["g"]:
                                                            for t_delay in p["t_delay"]:
                                                                params.append({
                                                                    "Ne": Ne,
                                                                    "tmin": tmin,
                                                                    "tmax": tmax,
                                                                    "Nt": p["Nt"],
                                                                    "T": T,
                                                                    "wd": wd,
                                                                    "s": p["s"],
                                                                    "m": m,
                                                                    "ef": ef,
                                                                    "g": g,
                                                                    "U": U,
                                                                    "A0": p["A0"],
                                                                    "tau": tau,
                                                                    "w":  w,
                                                                    "A0_pr": A0_pr,
                                                                    "tau_pr": tau_pr,
                                                                    "w_pr": w_pr,
                                                                    "t_delay": t_delay
                                                                })

print('parameters generated')
