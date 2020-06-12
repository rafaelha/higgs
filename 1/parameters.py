import numpy as np

params_ = [
    {
        "Ne": [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800],
        "tmin": [-4*(2*np.pi)],
        "tmax": [10*(2*np.pi)],
        "Nt": 200,
        "T": [0.02],
        "wd":  [10],
        "s": np.array([1,-1]),
        "m": [ np.array([1.0,1.2]) ],
        "ef": [ np.array([500,300]) ],
        "g": [ np.array([10, 10]) ],
        "U": [ np.array([[0.08,0.05], [0.05,0.18]]) ],
        "A0": 0.05,
        "tau": [2*np.pi],
        "w":  [1],
        "A0_pr": [0],
        "tau_pr": [0.1],
        "w_pr": [0],
        "t_delay": [0]
    }
]

params = []
for p in params_:
    for Ne in p["Ne"]:
        for tmin in p["tmin"]:
            for tmax in p["tmax"]:
                for T in p["T"]:
                    for wd in p["wd"]:
                        for m in p["m"]:
                            for ef in p["ef"]:
                                for g in p["g"]:
                                    for U in p["U"]:
                                        for tau in p["tau"]:
                                            for w in p["w"]:
                                                for A0_pr in p["A0_pr"]:
                                                    for tau_pr in p["tau_pr"]:
                                                        for w_pr in p["w_pr"]:
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