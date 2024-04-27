import numpy as np
import matplotlib.pyplot as plt

posOld = "E:\\Calculations\\Debugging Suite\\CompareWithPythonVersion3_samematerial_hourglass_it200_lam542_sym_filterOff_periodicFalse_beta0_epsilon_0.1_penaltyOff\\CoreStructure\\"
posNew = "E:\\Calculations\\Hourglass Python\\"

all_errors = [0]*200

for i in range(200):
    corestructureold = np.loadtxt(posOld + "CoreStructure" + str(i) + ".txt")
    corestructurenew = np.loadtxt(posNew + "CoreStructure" + str(i) + ".txt")

    for j in range(len(corestructurenew)):
        all_errors[i] += abs(corestructureold[j] - corestructurenew[j])

plt.figure(1)
plt.plot(all_errors)
#plt.yscale('log')
plt.title("L1 Norm")
plt.ylabel("Norm")
plt.xlabel("Iteration")
plt.show()