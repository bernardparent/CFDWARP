#!/usr/local/bin/soap

EoverN=20e-20; {Vm2}

EoverNtown=EoverN/1.0e-21;

{air}
Te_air=11606.0*(5.598*exp(0.000684*EoverNtown)-5.425*exp(-0.00281*EoverNtown));


{0.95 air + 0.05 C2H4}
Te_fuelair=11606.0*(5.446*exp(0.0007026*EoverNtown)-5.272*exp(-0.003088*EoverNtown));


printf("EoverN=%E  Vm2   Te_air=%E K    Te_fuelair=%E K\n",EoverN,Te_air,Te_fuelair);
