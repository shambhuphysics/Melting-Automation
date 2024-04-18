#########################################################
# Constants that you need to change  
# File Required #  Ueam_sol.dat Udft_sol.dat Pdft_sol.dat Ueam_liq.dat Udft_liq.dat Pdft_liq.dat
atomnum = 500
T = 867.41  
P=0.158  
K_sol = 45.49
vol_sol = 17.96
K_liq = 39.68
vol_liq = 18.88
delS_soliq = 0.0000967748736303918
std_error_delS = 0.00000125560396535062
###########################################################
import numpy as np
import pandas as pd
def calculate_free_energy(Ueam_data, Udft_data, Pdft_data, K, vol):
    # Calculate ΔU
    ediff = Udft_data - Ueam_data
    ediff.index = ediff.index + 1
    # Calculate total average U and its standard error
    avg_U = np.mean(ediff) / atomnum
    std_error_avg_U = np.std(ediff, ddof=1) / np.sqrt(len(ediff)) / atomnum

    # Calculate average squared difference and its standard error
    avg_squared_diff = np.mean((ediff - np.mean(ediff))**2) / atomnum
    std_error_avg_squared_diff = np.std((ediff - np.mean(ediff))**2, ddof=1) / np.sqrt(len(ediff))

    # Calculate beta * <ΔU^2> and its standard error
    beta_avgSq = avg_squared_diff / (2 * 8.61733e-05 * T)  #kb = 8.61733e-05  # Boltzmann constant in eV/K
    #std_error_beta_avgSq = std_error_avg_squared_diff / beta_avgSq * atomnum
    std_error_beta_avgSq = std_error_avg_squared_diff / (2 * atomnum*8.61733e-05*T)

    # Calculate total average U including its standard error
    total_avgU = avg_U - beta_avgSq
    #std_error_total_avgU = np.sqrt((std_error_avg_U)**2 + (avg_squared_diff)**2)
    std_error_total_avgU = np.sqrt((std_error_avg_U)**2 + (std_error_beta_avgSq)**2)
    # Calculate pressure-related values
    mean_pressure = np.mean(Pdft_data) / 10
    std_error_pres = np.std(Pdft_data, ddof=1) / np.sqrt(len(Pdft_data))
    delp_squ = (mean_pressure - P) ** 2

    # Calculate pressure contribution to solvation free energy and its standard error
    press_contrib = (vol * delp_squ) / (2 * K * 160.21766)
    std_error_press = ((vol * std_error_pres) / K) / 160.21766

    # Calculate solvation free energy and its standard error
    free_energy = total_avgU - press_contrib
    std_error_free_energy = np.sqrt((std_error_press)**2 + (std_error_total_avgU)**2)

    return free_energy, std_error_free_energy

def load_data(files):
    return [pd.read_csv(f, header=None, delim_whitespace=True).iloc[:, 1] for f in files]

# Solid data
solid_files = ['Ueam_sol.dat', 'Udft_sol.dat', 'Pdft_sol.dat']
Ueam_sol, Udft_sol, Pdft_sol = load_data(solid_files)

# Liquid data
liquid_files = ['Ueam_liq.dat', 'Udft_liq.dat', 'Pdft_liq.dat']
Ueam_liq, Udft_liq, Pdft_liq = load_data(liquid_files)

# Open the result file
with open('result.dat', 'w') as result_file:
    # Iterate over rows from 1 to the end of the .dat files with a step of 4
      Ueam_liq_data = pd.read_csv('Ueam_liq.dat', header=None, delim_whitespace=True)


      for num_rows in range(2, len(Ueam_liq_data),6):
        # Calculate solvation energy for solid
        G_sol, std_G_sol = calculate_free_energy(Ueam_sol[:num_rows], Udft_sol[:num_rows], Pdft_sol[:num_rows], K_sol, vol_sol)
        
        # Calculate solvation energy for liquid
        G_liq, std_G_liq = calculate_free_energy(Ueam_liq[:num_rows], Udft_liq[:num_rows], Pdft_liq[:num_rows], K_liq, vol_liq)

        # Calculate the difference in solvation free energy and its standard error
        DeltaG_soliq = G_liq - G_sol
        std_error_DeltaG_soliq = np.sqrt((std_G_sol)**2 + (std_G_liq)**2)

        # Calculate the change in temperature and its standard error
        DeltaT = DeltaG_soliq / delS_soliq
        std_error_DeltaT = DeltaT * np.sqrt((std_error_DeltaG_soliq / DeltaG_soliq)**2 + (std_error_delS / delS_soliq)**2)

        # Write the results to the result file
        print(f"{num_rows}, {DeltaT}, {std_error_DeltaT}", file=result_file)
        #print(f"{num_rows}, {DeltaG_soliq}, {std_error_DeltaG_soliq}", file=result_file)
        # Open and read the contents of the result file
with open('result.dat', 'r') as result_file:
    # Print each line of the file
    for line in result_file:
        print(line.strip())

import numpy as np
import matplotlib.pyplot as plt

# Read data from the result file
data = np.loadtxt('result.dat', delimiter=',', skiprows=0)

# Extract columns for plotting
num_rows = data[:, 0]
DeltaT = data[:, 1]
std_error_DeltaT = data[:, 2]

# Plot
plt.errorbar(num_rows, DeltaT, yerr=std_error_DeltaT, fmt='o', color='blue', ecolor='red', capsize=5)

# Add lines between the points
plt.plot(num_rows, DeltaT, linestyle='-', color='blue')

plt.xlabel('Number of Rows')
plt.ylabel('DeltaT')
plt.title('DeltaT vs Number of Rows')
plt.grid(True)

# Save the plot in PNG format
plt.savefig('plot_with_lines.png')

# Save the plot in JPG format
plt.savefig('plot_with_lines.jpg')

# Save the plot in PDF format
plt.savefig('plot_with_lines.pdf')

plt.show()

