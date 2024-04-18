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
from datetime import datetime

# Get the current date
current_date = datetime.now().strftime("%Y-%m-%d")

print("   ********* Welcome to Melting Calculator ******** ")
print(f"           Date of Calculation: {current_date}")
print("   ********* Results for Solid ******** ")

# Solid data
solid_files = ['Ueam_sol.dat', 'Udft_sol.dat', 'Pdft_sol.dat']
# Liquid data
liquid_files = ['Ueam_liq.dat', 'Udft_liq.dat', 'Pdft_liq.dat']


def calculate_Free_Energy(Ueam_data, Udft_data, Pdft_data, K, vol):
    # Calculate ΔU
    ediff = Udft_data - Ueam_data 
    ediff.index = ediff.index + 1  # Adjust index to start from 1

    #print('ediff', ediff)  # good

    # Calculate total average U and its standard error
    avg_U = np.mean(ediff) / atomnum
    std_error_avg_U = np.std(ediff, ddof=1) / np.sqrt(len(ediff)) / atomnum
    print('  <ΔU/N> =', avg_U,'±',std_error_avg_U,'eV/atom')  # good
    avg_squared_diff = np.mean((ediff - np.mean(ediff))**2) / atomnum
    std_error_avg_squared_diff = np.std((ediff - np.mean(ediff))**2, ddof=1) / np.sqrt(len(ediff))
    # Calculate beta * <ΔU^2> and its standard error
    beta_avgSq = avg_squared_diff / (2 * 8.61733e-05 * T)
    std_error_beta_avgSq = std_error_avg_squared_diff / (2 * atomnum*8.61733e-05*T)

    # Calculate total average U including its standard error
    total_avgU = avg_U - beta_avgSq
    std_error_total_avgU = np.sqrt((std_error_avg_U)**2 + (std_error_beta_avgSq)**2) #std_error_beta_avgSq
    print('  <ΔF/N> =', total_avgU,'±',std_error_total_avgU,'eV/atom')
    
    
    # Calculate pressure-related values
    mean_pressure = np.mean(Pdft_data) / 10
    std_error_pres = np.std(Pdft_data, ddof=1) / np.sqrt(len(Pdft_data))
    delp_squ = (mean_pressure -P) ** 2
    # Calculate pressure contribution to free energy and its standard error
    press_contrib = (vol * delp_squ) / (2 * K * 160.21766)
    std_error_press = ((vol * std_error_pres) / K) / 160.21766
    # Calculate free energy and its standard error
    Free_Energy = total_avgU - press_contrib
    std_error_Free_Energy = np.sqrt((std_error_press)**2 + (std_error_total_avgU)**2)
    print('  <ΔG/N> =', Free_Energy,' ±',std_error_Free_Energy,'eV/atom')
    return Free_Energy, std_error_Free_Energy



# Load solid data
Ueam_sol, Udft_sol, Pdft_sol = [pd.read_csv(f, header=None, delim_whitespace=True).iloc[:, 1] for f in solid_files]

# Load liquid data
Ueam_liq, Udft_liq, Pdft_liq = [pd.read_csv(f, header=None, delim_whitespace=True).iloc[:, 0] + 1 for f in liquid_files]
#Ueam_liq, Udft_liq, Pdft_liq = [pd.read_csv(f, header=None, delim_whitespace=True).iloc[:, 0] + 1 for f in liquid_files]

# Calculate energy for solid
G_sol, std_G_sol = calculate_Free_Energy(Ueam_sol, Udft_sol, Pdft_sol, K_sol, vol_sol)

print("   ********* Results for Liquid ******** ")
#print("Standard Error of Free Energy (Solid):", std_G_sol)

# Load liquid data
Ueam_liq, Udft_liq, Pdft_liq = [pd.read_csv(f, header=None, delim_whitespace=True).iloc[:, 1] for f in liquid_files]

# Calculate energy for liquid
G_liq, std_G_liq = calculate_Free_Energy(Ueam_liq, Udft_liq, Pdft_liq, K_liq, vol_liq)

#print("Energy (Liquid):", G_liq)
#print("Standard Error of Free Energy (Liquid):", std_G_liq)
print("   ********* Result for Liquid - Solid ******** ")
DeltaG_soliq = G_liq -G_sol
std_error_DeltaG_soliq=np.sqrt((std_G_sol)**2+(std_G_liq)**2)
print ('    ΔG_ls =', DeltaG_soliq, '±', std_error_DeltaG_soliq,'eV/atom')

DeltaT= DeltaG_soliq/delS_soliq
std_error_DeltaT=DeltaT*np.sqrt((std_error_DeltaG_soliq/DeltaG_soliq)**2+(std_error_delS/delS_soliq)**2)
print ('    ΔT =', DeltaT, '±', std_error_DeltaT,'K')
print(f"Hence, The ab-initio Tm at {P} GPa = {DeltaT+T} ± {std_error_DeltaT}")
