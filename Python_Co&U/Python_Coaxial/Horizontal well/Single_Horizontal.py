import numpy as np

def twat_lateral_time_evaluation(Tr, Ti, Kr, L, mflow, Cpf, alpha_r, t, rw):
    Tout = Tr - ((Tr - Ti)/(np.exp((2 * np.pi * Kr * L)/(mflow * Cpf * np.log((np.pi * (alpha_r * t)**0.5)/(2 * rw))))))
    return Tout 


if __name__ == "__main__":
    Twat_out =  twat_lateral_time_evaluation(77.8, 24, 4.64, 1900, 2.88, 4180, 1.57 * 10**(-6), 31536000, 0.079)
    print(Twat_out)