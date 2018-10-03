#!/bin/csh


# Dominique Marion, Diego F. Gauto, Isabel Ayala, Audrey Hessel, Paul Schanda 
# Univ. Grenoble Alpes, CEA CNRS, Institute for Structural Biology (IBS) 71 avenue des martyrs,
# 38044 Grenoble (France); E-mail: dominique.marion@ibs.fr, paul.schanda@ibs.fr

# Simulation of the t1rho values at 600MHz for Ubiquitin on Agilent spectrometer


# Chemical shift difference    -> 16 values
# Dipolar angular fluctuation  -> 13 values
# Population                   -> 8 values
# Exchange rate                -> 21 values
# Rf fields                    -> 28 values
rm Summary


setenv J_coupling 90.0
# Dipolar coupling with S=0.98
setenv D_coupling 22495
setenv D_alpha 0.0
setenv D_beta 103.0
setenv D_gamma 0.0

setenv N_shift 0
setenv H_shift 0

# CSA for 600 MHz with S=0.98
setenv N_CSA_delta 6600
setenv N_CSA_eta 0  
setenv N_CSA_alpha 0  
setenv N_CSA_beta 83  
setenv N_CSA_gamma 0  

setenv H_CSA_delta 480
setenv H_CSA_eta 0  
setenv H_CSA_alpha 0  
setenv H_CSA_beta 0  
setenv H_CSA_gamma 0  

setenv Powder_quality 2
setenv Time_steps 100
setenv Spin_speed 35000
setenv Sampling_rate 4
setenv Relax_delay 0.15
setenv Random_field 0

setenv Orientation 0
setenv Coherence -1
setenv Decoupling 1
setenv Hstate 0

echo -n "# Spinning_speed " >> Summary
printenv Spin_speed >> Summary
echo -n "# Relaxation_delay " >> Summary
printenv  Relax_delay >> Summary
echo "# [Input] Population  Exch_rate RF_strength  Chem_shift_diff  Dip_tensor_angle" >> Summary
echo "# [Output] MonoExp BiExp[1] Biexp[2] Ampl[1] Ampl[2] Filename" >> Summary
echo "#     [Input]     <- | ->    [Output] ">> Summary

@ counter = 0

# Population list: 2 5 10 15 20 25 30 50

@ population = 5
setenv Spin_Sys_Name 2spin_5perc.sys 


#--- First simulation for 600 MHz ---

#--- Chemical exchange --------------

@ chem_dif = 1000
@ dip_angle = 0

foreach ch (0 4 8 12 16 20 25 30 35 40 45 50 60 70 80 100)

foreach exch (0 125 500 1125 2000 3125 4500 6125 8000 10125 12500 15125 18000 21125 24500 28125 32000 36125 40500 45125 50000)



# rf field values for Agilent experiments
# foreach i (9 10 12 14 16 18 20 22 24 32 40 48 60 80 88 96 104 112 120 124 128 130 132 )
foreach i (7 8 9 10 11 12 14 16 18 20 22 24 32 40 48 60 80 88 96 104 112 120 124 128 129 130 131 132 )


@ rf_field = 250 * $i 
@ kex =  $exch
# 600 MHZ spectrometer
@ chem_dif = 6 * $ch
@ counter = $counter + 1


set filename = Simul_310118_Nx_${counter}_${population}
echo -n ${population} ${kex} ${rf_field} ${chem_dif} ${dip_angle} "    " >> Summary 

./R1rho_Components_pop  ${kex} ${rf_field} ${chem_dif} ${dip_angle} ${filename}
python Fit_Mono_Exp.py ${filename} "Nx kex "${kex} 


end
end
end

@ chem_dif = 0
@ dip_angle = 20

#--- Dipolar fluctuation --------------

foreach an (0 2 4 6 8 10 12 14 16 20 24 30 36)

foreach exch (0 125 500 1125 2000 3125 4500 6125 8000 10125 12500 15125 18000 21125 24500 28125 32000 36125 40500 45125 50000)

foreach i (7 8 9 10 11 12 14 16 18 20 22 24 32 40 48 60 80 88 96 104 112 120 124 128 129 130 131 132 )


@ rf_field = 250 * $i 
@ kex =  $exch
@ dip_angle = $an
@ counter = $counter + 1

set filename = Simul_310118_Nx_${counter}_${population}
echo -n ${population} ${kex} ${rf_field} ${chem_dif} ${dip_angle} "    " >> Summary 

./R1rho_Components_pop  ${kex} ${rf_field} ${chem_dif} ${dip_angle} ${filename}
python Fit_Mono_Exp.py ${filename} "Nx kex "${kex} 


end
end
end


cp Summary Summary_env_Simul_310118_${population}

