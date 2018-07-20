#!/bin/sh
echo "Removing proprietary directory .././model/chem/_Air_Plasma_8s_28r_Macheret"
rm -rf .././model/chem/_Air_Plasma_8s_28r_Macheret 
cp -a _proprietary .././model/chem/_Air_Plasma_8s_28r_Macheret 
( printf "current Air_Plasma_8s_28r_Macheret\nEND" ) > .././model/chem/_Air_Plasma_8s_28r_Macheret/.config 
echo "Removing proprietary directory .././model/chem/_H2_Air_Plasma_14s_47r_Jachimowsky-Macheret"
rm -rf .././model/chem/_H2_Air_Plasma_14s_47r_Jachimowsky-Macheret 
cp -a _proprietary .././model/chem/_H2_Air_Plasma_14s_47r_Jachimowsky-Macheret 
( printf "current H2_Air_Plasma_14s_47r_Jachimowsky-Macheret\nEND" ) > .././model/chem/_H2_Air_Plasma_14s_47r_Jachimowsky-Macheret/.config 
echo "Removing proprietary directory .././model/fluid/_Navier-Stokes_plasma"
rm -rf .././model/fluid/_Navier-Stokes_plasma 
cp -a _proprietary .././model/fluid/_Navier-Stokes_plasma 
( printf "current Navier-Stokes_plasma\nEND" ) > .././model/fluid/_Navier-Stokes_plasma/.config 
echo "Removing proprietary directory .././model/fluid/_Favre-Reynolds_plasma"
rm -rf .././model/fluid/_Favre-Reynolds_plasma 
cp -a _proprietary .././model/fluid/_Favre-Reynolds_plasma 
( printf "current Favre-Reynolds_plasma\nEND" ) > .././model/fluid/_Favre-Reynolds_plasma/.config 
echo "Removing proprietary directory .././model/fluid/_drift-diffusion"
rm -rf .././model/fluid/_drift-diffusion 
cp -a _proprietary .././model/fluid/_drift-diffusion 
( printf "current drift-diffusion\nEND" ) > .././model/fluid/_drift-diffusion/.config 
echo "Removing proprietary directory .././model/emfield/_Ohm"
rm -rf .././model/emfield/_Ohm 
cp -a _proprietary .././model/emfield/_Ohm 
( printf "current Ohm\nEND" ) > .././model/emfield/_Ohm/.config 
echo "Removing proprietary directory .././model/emfield/_Ohm_generalized"
rm -rf .././model/emfield/_Ohm_generalized 
cp -a _proprietary .././model/emfield/_Ohm_generalized 
( printf "current Ohm_generalized\nEND" ) > .././model/emfield/_Ohm_generalized/.config 
echo "Removing proprietary directory .././cycle/ts/_block_DDADI"
rm -rf .././cycle/ts/_block_DDADI 
cp -a _proprietary .././cycle/ts/_block_DDADI 
( printf "current block_DDADI\nEND" ) > .././cycle/ts/_block_DDADI/.config 
echo "Removing proprietary directory .././cycle/ts/_block_IMAF"
rm -rf .././cycle/ts/_block_IMAF 
cp -a _proprietary .././cycle/ts/_block_IMAF 
( printf "current block_IMAF\nEND" ) > .././cycle/ts/_block_IMAF/.config 
echo "Removing proprietary directory .././cycle/tsemf/_AF"
rm -rf .././cycle/tsemf/_AF 
cp -a _proprietary .././cycle/tsemf/_AF 
( printf "current AF\nEND" ) > .././cycle/tsemf/_AF/.config 
