// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2023 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>



/* taken from  Zettervall, N., Fureby, C., and Nilsson, E. J., 
  “Small skeletal kinetic reaction mechanism for ethylene–air combustion,” 
   Energy and Fuels, Vol. 31, No. 12, 2017, pp. 14138–14149. */

const static bool REACTION[67] = {
  TRUE, /* reaction[0] */
  TRUE, /* reaction[1] */
  TRUE, /* reaction[2] */
  TRUE, /* reaction[3] */
  TRUE, /* reaction[4] */
  TRUE, /* reaction[5] */
  TRUE, /* reaction[6] */
  TRUE, /* reaction[7] */
  TRUE, /* reaction[8] */
  TRUE, /* reaction[9] */
  TRUE, /* reaction[10] */
  TRUE, /* reaction[11] */
  TRUE, /* reaction[12] */
  TRUE, /* reaction[13] */
  TRUE, /* reaction[14] */
  TRUE, /* reaction[15] */
  TRUE, /* reaction[16] */
  TRUE, /* reaction[17] */
  TRUE, /* reaction[18] */
  TRUE, /* reaction[19] */
  TRUE, /* reaction[20] */
  TRUE, /* reaction[21] */
  TRUE, /* reaction[22] */
  TRUE, /* reaction[23] */
  TRUE, /* reaction[24] */
  TRUE, /* reaction[25] */
  TRUE, /* reaction[26] */
  TRUE, /* reaction[27] */
  TRUE, /* reaction[28] */
  TRUE, /* reaction[29] */
  TRUE, /* reaction[30] */
  TRUE, /* reaction[31] */
  TRUE, /* reaction[32] */
  TRUE, /* reaction[33] */
  TRUE, /* reaction[34] */
  TRUE, /* reaction[35] */
  TRUE, /* reaction[36] */
  TRUE, /* reaction[37] */
  TRUE, /* reaction[38] */
  TRUE, /* reaction[39] */
  TRUE, /* reaction[40] */
  TRUE, /* reaction[41] */
  TRUE, /* reaction[42] */
  TRUE, /* reaction[43] */
  TRUE, /* reaction[44] */
  TRUE, /* reaction[45] */
  TRUE, /* reaction[46] */
  TRUE, /* reaction[47] */
  TRUE, /* reaction[48] */
  TRUE, /* reaction[49] */
  TRUE, /* reaction[50] */
  TRUE, /* reaction[51] */
  TRUE, /* reaction[52] */
  TRUE, /* reaction[53] */
  TRUE, /* reaction[54] */
  TRUE, /* reaction[55] */
  TRUE, /* reaction[56] */
  TRUE, /* reaction[57] */
  TRUE, /* reaction[58] */
  TRUE, /* reaction[59] */
  TRUE, /* reaction[60] */
  TRUE, /* reaction[61] */
  TRUE, /* reaction[62] */
  TRUE, /* reaction[63] */
  TRUE, /* reaction[64] */
  TRUE, /* reaction[65] */
  TRUE, /* reaction[66] */
};

void find_W_Zettervall2017b ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  double eff;
  long k, specM;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM (k) * 1.0e-6;    /* mole/cm^3 */
    W[k] = 0.0;
  }

  if (REACTION[1])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specC2H4, specH, specM, specC2H5, specM, 4.17e+10, 0.0, 11030.0, T, X, W);
    }
  }

  if (REACTION[2])
    add_to_W_fw_2r2p(specC2H5, specH, specCH3, specCH3, 1.5e+13, 0.0, 0.0, T, X, W);

  if (REACTION[3])
    add_to_W_fw_1r2p(specC2H4, specC2H3, specH, 1.0e+12, 0.0, 48200.0, T, X, W);

  if (REACTION[4])
    add_to_W_fw_2r1p(specC2H4, specH, specC2H5, 1.8e+12, 1.0, 9400.0, T, X, W);

  if (REACTION[5])
    add_to_W_fw_2r2p(specC2H4, specO, specHCO, specCH3, 3.31e+12, 0.0, 1130.0, T, X, W);

  if (REACTION[6])
    add_to_W_fw_2r2p(specHCO, specCH3, specC2H4, specO, 1.58e+11, 0.0, 31180.0, T, X, W);

  if (REACTION[7])
    add_to_W_fw_2r2p(specC2H4, specOH, specC2H3, specH2O, 4.79e+12, 0.0, 1230.0, T, X, W);

  if (REACTION[8])
    add_to_W_fw_2r2p(specC2H3, specH2O, specC2H4, specOH, 1.2e+12, 0.0, 14000.0, T, X, W);

  if (REACTION[9])
    add_to_W_fw_2r2p(specC2H4, specCH3, specC2H3, specCH4, 1.0e+13, 0.0, 13000.0, T, X, W);

  if (REACTION[10])
    add_to_W_fw_2r2p(specC2H3, specCH4, specC2H4, specCH3, 3.02e+13, 0.0, 12580.0, T, X, W);

  if (REACTION[11])
    add_to_W_fw_2r2p(specC2H4, specH, specC2H3, specH2, 5.0e+12, 0.0, 5000.0, T, X, W);

  if (REACTION[12])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specC2H2, specH, specM, specC2H3, specM, 1.23e+11, 1.0, 10360.0, T, X, W);
    }
  }

  if (REACTION[13])
    add_to_W_fw_2r2p(specC2H3, specH, specC2H2, specH2, 2.0e+13, 0.0, 2500.0, T, X, W);

  if (REACTION[14])
    add_to_W_fw_2r2p(specC2H3, specOH, specC2H2, specH2O, 7.0e+14, 0.0, 1100.0, T, X, W);

  if (REACTION[15])
    add_to_W_fw_1r2p(specC2H3, specC2H, specH2, 1.8e+12, 0.0, 44800.0, T, X, W);

  if (REACTION[16])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specC2H, specH, specM, specC2H2, specM, 1.1e+09, 1.0, 770.0, T, X, W);
    }
  }

  if (REACTION[17])
    add_to_W_fw_2r2p(specC2H2, specH, specC2H, specH2, 2.0e+14, 0.0, 19000.0, T, X, W);

  if (REACTION[18])
    add_to_W_fw_2r2p(specC2H2, specOH, specC2H, specH2O, 8.0e+12, 0.0, 5000.0, T, X, W);

  if (REACTION[19])
    add_to_W_fw_2r2p(specC2H, specH2O, specC2H2, specOH, 5.37e+12, 0.0, 16360.0, T, X, W);

  if (REACTION[20])
    add_to_W_fw_2r2p(specC2H2, specO, specC2H, specOH, 3.24e+15, 0.6, 12000.0, T, X, W);

  if (REACTION[21])
    add_to_W_fw_2r2p(specC2H, specOH, specC2H2, specO, 2.95e+14, 0.6, 910.0, T, X, W);

  if (REACTION[22])
    add_to_W_fw_2r2p(specC2H, specO2, specHCO, specCO, 1.0e+13, 0.0, 7000.0, T, X, W);

  if (REACTION[23])
    add_to_W_fw_2r2p(specHCO, specCO, specC2H, specO2, 8.51e+12, 0.0, 138400.0, T, X, W);

  if (REACTION[24])
    add_to_W_fw_2r2p(specC2H, specOH, specHCO, specCH, 8.0e+12, 0.0, 5000.0, T, X, W);

  if (REACTION[25])
    add_to_W_fw_2r2p(specH, specO2, specOH, specO, 5.0E+14, 0.0, 16800.0, T, X, W);

  if (REACTION[26])
    add_to_W_fw_2r2p(specOH, specO, specH, specO2, 1.20E+13, 0.0, 690.0, T, X, W);

  if (REACTION[27])
    add_to_W_fw_2r2p(specO, specH2, specOH, specH, 1.80000E+10, 1.0, 8826.0, T, X, W);

  if (REACTION[28])
    add_to_W_fw_2r2p(specOH, specH, specO, specH2, 8.00000E+09, 1.0, 6760.0, T, X, W);

  if (REACTION[29])
    add_to_W_fw_2r2p(specH2, specOH, specH2O, specH, 1.17000E+09, 1.3, 3626.0, T, X, W);

  if (REACTION[30])
    add_to_W_fw_2r2p(specH2O, specH, specH2, specOH, 6.0E+09, 1.3, 18588.0, T, X, W);

  if (REACTION[31])
    add_to_W_fw_2r2p(specOH, specOH, specO, specH2O, 6.00000E+08, 1.3, 0.0, T, X, W);

  if (REACTION[32])
    add_to_W_fw_2r2p(specO, specH2O, specOH, specOH, 4.0E+09, 1.3, 17029.0, T, X, W);

  if (REACTION[33])  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        case specO2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_W_fw_3r2p(specH, specO2, specM, specHO2, specM, eff*2.5E+18, -0.8, 0.0, T, X, W);
    }
  }

  if (REACTION[34])
    add_to_W_fw_2r2p(specH, specHO2, specOH, specOH, 1.50000E+14, 0.0, 1004.0, T, X, W);

  if (REACTION[35])
    add_to_W_fw_2r2p(specH, specHO2, specH2, specO2, 2.50000E+13, 0.0, 700.0, T, X, W);

  if (REACTION[36])
    add_to_W_fw_2r2p(specOH, specHO2, specH2O, specO2, 2.00000E+13, 0.0, 1000.0, T, X, W);

  if (REACTION[37])
    add_to_W_fw_2r2p(specCO, specOH, specCO2, specH, 1.51000E+07, 1.3, -758.0, T, X, W);

  if (REACTION[38])
    add_to_W_fw_2r2p(specCO2, specH, specCO, specOH, 1.57000E+09, 1.3, 19800.0, T, X, W);

  if (REACTION[39] && gl->model.chem.LINDEMANNREACTIONS)  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specC2H4: eff = 3.0; break;
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        case specO2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_W_fw_2r3p_Lindemann(specCH4, specM, specCH3, specH, specM, eff*6.30000E+14, 0.0, 104000.0, eff*1.00000E+17, 0.0, 86000.0, T, X, W);
    }
  }

  if (REACTION[40] && gl->model.chem.LINDEMANNREACTIONS)  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specC2H4: eff = 3.0; break;
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_W_fw_3r2p_Lindemann(specCH3, specH, specM, specCH4, specM, eff*5.20000E+12, 0.0, -1310.0, eff*8.25000E+14, 0.0, -19310.0, max(gl->model.chem.TMIN_LINDEMANN,T), X, W);
    }
  }

  if (REACTION[41])
    add_to_W_fw_2r2p(specCH4, specH, specCH3, specH2, 2.20000E+04, 3.0, 8750.0, T, X, W);

  if (REACTION[42])
    add_to_W_fw_2r2p(specCH3, specH2, specCH4, specH, 9.57000E+02, 3.0, 8750.0, T, X, W);

  if (REACTION[43])
    add_to_W_fw_2r2p(specCH4, specOH, specCH3, specH2O, 1.60000E+06, 2.1, 2460.0, T, X, W);

  if (REACTION[44])
    add_to_W_fw_2r2p(specCH3, specH2O, specCH4, specOH, 3.02000E+05, 2.1, 17422.0, T, X, W);

  if (REACTION[45])
    add_to_W_fw_2r2p(specCH3, specO, specCH2O, specH, 6.80000E+13, 0.0, 0.0, T, X, W);

  if (REACTION[46])
    add_to_W_fw_2r2p(specCH2O, specH, specHCO, specH2, 5.0E+13, 0.0, 3991.0, T, X, W);

  if (REACTION[47])
    add_to_W_fw_2r2p(specCH2O, specOH, specHCO, specH2O, 1.4E+14, 0.0, 1100.0, T, X, W);

  if (REACTION[48])
    add_to_W_fw_2r2p(specHCO, specH, specCO, specH2, 4.00000E+13, 0.0, 0.0, T, X, W);

  if (REACTION[49])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_2r3p(specHCO, specM, specCO, specH, specM, 1.60000E+14, 0.0, 14700.0, T, X, W);
    }
  }

  if (REACTION[50])
    add_to_W_fw_2r2p(specCH3, specO2, specCH3O, specO, 3.0E+13, 0.0, 25652.0, T, X, W);

  if (REACTION[51])
    add_to_W_fw_2r2p(specCH3O, specH, specCH2O, specH2, 2.00000E+13, 0.0, 0.0, T, X, W);

  if (REACTION[52])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_2r3p(specCH3O, specM, specCH2O, specH, specM, 2.40000E+13, 0.0, 28812.0, T, X, W);
    }
  }

  if (REACTION[53])
    add_to_W_fw_2r2p(specHO2, specHO2, specH2O2, specO2, 2.0E+14, 0.0, 0.0, T, X, W);

  if (REACTION[54])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_2r3p(specH2O2, specM, specOH, specOH, specM, 1.30000E+17, 0.0, 45500.0, T, X, W);
    }
  }

  if (REACTION[55])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specOH, specOH, specM, specH2O2, specM, 9.86000E+14, 0.0, -5070.0, T, X, W);
    }
  }

  if (REACTION[56])
    add_to_W_fw_2r2p(specH2O2, specOH, specH2O, specHO2, 1.00000E+13, 0.0, 1800.0, T, X, W);

  if (REACTION[57])
    add_to_W_fw_2r2p(specH2O, specHO2, specH2O2, specOH, 2.86000E+13, 0.0, 32790.0, T, X, W);

  if (REACTION[58])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specOH, specH, specM, specH2O, specM, 2.20000E+22, -2.0, 0.0, T, X, W);
    }
  }

  if (REACTION[59])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_W_fw_3r2p(specH, specH, specM, specH2, specM, 1.80000E+18, -1.0, 0.0, T, X, W);
    }
  }

  if (REACTION[60])
    add_to_W_fw_2r2p(specCH3, specOH, specCH2, specH2O, 7.6E+06, 2.0, 5000.0, T, X, W);

  if (REACTION[61])
    add_to_W_fw_2r2p(specCH2, specO, specCO, specH2, 3.0E+13, 0.0, 0.0, T, X, W);

  if (REACTION[62])
    add_to_W_fw_2r2p(specCH2, specOH, specCH, specH2O, 1.13E+07, 2.0, 3000.0, T, X, W);

  if (REACTION[63])
    add_to_W_fw_2r2p(specCH, specO, specCO, specH, 5.7E+13, 0.0, 0.0, T, X, W);

  if (REACTION[64])
    add_to_W_fw_2r2p(specCH, specOH, specHCO, specH, 3.0E+13, 0.0, 0.0, T, X, W);

  if (REACTION[65])
    add_to_W_fw_2r2p(specCH, specO2, specHCO, specO, 3.3E+13, 0.0, 0.0, T, X, W);

  if (REACTION[66])  
    add_to_W_fw_2r2p ( specCH, specCO2, specHCO, specCO,  8.40e+13, 0.0, 200.0, T, X, W );   

}

void find_dW_dx_Zettervall2017b ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long s, k, specM;
  double eff;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM (k) * 1.0e-6;    /* mole/cm^3 */
  }

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++) {
      dWdrhok[s][k] = 0.0;
    }
  }

  if (REACTION[1])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specC2H4, specH, specM, specC2H5, specM, 4.17e+10, 0.0, 11030.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[2])
    add_to_dW_fw_2r2p(specC2H5, specH, specCH3, specCH3, 1.5e+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[3])
    add_to_dW_fw_1r2p(specC2H4, specC2H3, specH, 1.0e+12, 0.0, 48200.0, T, X, dWdT, dWdrhok);

  if (REACTION[4])
    add_to_dW_fw_2r1p(specC2H4, specH, specC2H5, 1.8e+12, 1.0, 9400.0, T, X, dWdT, dWdrhok);

  if (REACTION[5])
    add_to_dW_fw_2r2p(specC2H4, specO, specHCO, specCH3, 3.31e+12, 0.0, 1130.0, T, X, dWdT, dWdrhok);

  if (REACTION[6])
    add_to_dW_fw_2r2p(specHCO, specCH3, specC2H4, specO, 1.58e+11, 0.0, 31180.0, T, X, dWdT, dWdrhok);

  if (REACTION[7])
    add_to_dW_fw_2r2p(specC2H4, specOH, specC2H3, specH2O, 4.79e+12, 0.0, 1230.0, T, X, dWdT, dWdrhok);

  if (REACTION[8])
    add_to_dW_fw_2r2p(specC2H3, specH2O, specC2H4, specOH, 1.2e+12, 0.0, 14000.0, T, X, dWdT, dWdrhok);

  if (REACTION[9])
    add_to_dW_fw_2r2p(specC2H4, specCH3, specC2H3, specCH4, 1.0e+13, 0.0, 13000.0, T, X, dWdT, dWdrhok);

  if (REACTION[10])
    add_to_dW_fw_2r2p(specC2H3, specCH4, specC2H4, specCH3, 3.02e+13, 0.0, 12580.0, T, X, dWdT, dWdrhok);

  if (REACTION[11])
    add_to_dW_fw_2r2p(specC2H4, specH, specC2H3, specH2, 5.0e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok);

  if (REACTION[12])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specC2H2, specH, specM, specC2H3, specM, 1.23e+11, 1.0, 10360.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[13])
    add_to_dW_fw_2r2p(specC2H3, specH, specC2H2, specH2, 2.0e+13, 0.0, 2500.0, T, X, dWdT, dWdrhok);

  if (REACTION[14])
    add_to_dW_fw_2r2p(specC2H3, specOH, specC2H2, specH2O, 7.0e+14, 0.0, 1100.0, T, X, dWdT, dWdrhok);

  if (REACTION[15])
    add_to_dW_fw_1r2p(specC2H3, specC2H, specH2, 1.8e+12, 0.0, 44800.0, T, X, dWdT, dWdrhok);

  if (REACTION[16])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specC2H, specH, specM, specC2H2, specM, 1.1e+09, 1.0, 770.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[17])
    add_to_dW_fw_2r2p(specC2H2, specH, specC2H, specH2, 2.0e+14, 0.0, 19000.0, T, X, dWdT, dWdrhok);

  if (REACTION[18])
    add_to_dW_fw_2r2p(specC2H2, specOH, specC2H, specH2O, 8.0e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok);

  if (REACTION[19])
    add_to_dW_fw_2r2p(specC2H, specH2O, specC2H2, specOH, 5.37e+12, 0.0, 16360.0, T, X, dWdT, dWdrhok);

  if (REACTION[20])
    add_to_dW_fw_2r2p(specC2H2, specO, specC2H, specOH, 3.24e+15, 0.6, 12000.0, T, X, dWdT, dWdrhok);

  if (REACTION[21])
    add_to_dW_fw_2r2p(specC2H, specOH, specC2H2, specO, 2.95e+14, 0.6, 910.0, T, X, dWdT, dWdrhok);

  if (REACTION[22])
    add_to_dW_fw_2r2p(specC2H, specO2, specHCO, specCO, 1.0e+13, 0.0, 7000.0, T, X, dWdT, dWdrhok);

  if (REACTION[23])
    add_to_dW_fw_2r2p(specHCO, specCO, specC2H, specO2, 8.51e+12, 0.0, 138400.0, T, X, dWdT, dWdrhok);

  if (REACTION[24])
    add_to_dW_fw_2r2p(specC2H, specOH, specHCO, specCH, 8.0e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok);

  if (REACTION[25])
    add_to_dW_fw_2r2p(specH, specO2, specOH, specO, 5.0E+14, 0.0, 16800.0, T, X, dWdT, dWdrhok);

  if (REACTION[26])
    add_to_dW_fw_2r2p(specOH, specO, specH, specO2, 1.20E+13, 0.0, 690.0, T, X, dWdT, dWdrhok);

  if (REACTION[27])
    add_to_dW_fw_2r2p(specO, specH2, specOH, specH, 1.80000E+10, 1.0, 8826.0, T, X, dWdT, dWdrhok);

  if (REACTION[28])
    add_to_dW_fw_2r2p(specOH, specH, specO, specH2, 8.00000E+09, 1.0, 6760.0, T, X, dWdT, dWdrhok);

  if (REACTION[29])
    add_to_dW_fw_2r2p(specH2, specOH, specH2O, specH, 1.17000E+09, 1.3, 3626.0, T, X, dWdT, dWdrhok);

  if (REACTION[30])
    add_to_dW_fw_2r2p(specH2O, specH, specH2, specOH, 6.0E+09, 1.3, 18588.0, T, X, dWdT, dWdrhok);

  if (REACTION[31])
    add_to_dW_fw_2r2p(specOH, specOH, specO, specH2O, 6.00000E+08, 1.3, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[32])
    add_to_dW_fw_2r2p(specO, specH2O, specOH, specOH, 4.0E+09, 1.3, 17029.0, T, X, dWdT, dWdrhok);

  if (REACTION[33])  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        case specO2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_dW_fw_3r2p(specH, specO2, specM, specHO2, specM, eff*2.5E+18, -0.8, 0.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[34])
    add_to_dW_fw_2r2p(specH, specHO2, specOH, specOH, 1.50000E+14, 0.0, 1004.0, T, X, dWdT, dWdrhok);

  if (REACTION[35])
    add_to_dW_fw_2r2p(specH, specHO2, specH2, specO2, 2.50000E+13, 0.0, 700.0, T, X, dWdT, dWdrhok);

  if (REACTION[36])
    add_to_dW_fw_2r2p(specOH, specHO2, specH2O, specO2, 2.00000E+13, 0.0, 1000.0, T, X, dWdT, dWdrhok);

  if (REACTION[37])
    add_to_dW_fw_2r2p(specCO, specOH, specCO2, specH, 1.51000E+07, 1.3, -758.0, T, X, dWdT, dWdrhok);

  if (REACTION[38])
    add_to_dW_fw_2r2p(specCO2, specH, specCO, specOH, 1.57000E+09, 1.3, 19800.0, T, X, dWdT, dWdrhok);

  if (REACTION[39] && gl->model.chem.LINDEMANNREACTIONS)  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specC2H4: eff = 3.0; break;
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        case specO2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_dW_fw_2r3p_Lindemann(specCH4, specM, specCH3, specH, specM, eff*6.30000E+14, 0.0, 104000.0, eff*1.00000E+17, 0.0, 86000.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[40] && gl->model.chem.LINDEMANNREACTIONS)  {
    for (specM = 0; specM < ns; specM++) {
      switch (specM) {
        case specC2H4: eff = 3.0; break;
        case specCH4: eff = 6.5; break;
        case specCO: eff = 0.75; break;
        case specCO2: eff = 1.5; break;
        case specH2: eff = 1.0; break;
        case specH2O: eff = 6.5; break;
        case specN2: eff = 0.4; break;
        default: eff = 1.0;
      }
      add_to_dW_fw_3r2p_Lindemann(specCH3, specH, specM, specCH4, specM, eff*5.20000E+12, 0.0, -1310.0, eff*8.25000E+14, 0.0, -19310.0, max(gl->model.chem.TMIN_LINDEMANN,T), X, dWdT, dWdrhok);
    }
  }

  if (REACTION[41])
    add_to_dW_fw_2r2p(specCH4, specH, specCH3, specH2, 2.20000E+04, 3.0, 8750.0, T, X, dWdT, dWdrhok);

  if (REACTION[42])
    add_to_dW_fw_2r2p(specCH3, specH2, specCH4, specH, 9.57000E+02, 3.0, 8750.0, T, X, dWdT, dWdrhok);

  if (REACTION[43])
    add_to_dW_fw_2r2p(specCH4, specOH, specCH3, specH2O, 1.60000E+06, 2.1, 2460.0, T, X, dWdT, dWdrhok);

  if (REACTION[44])
    add_to_dW_fw_2r2p(specCH3, specH2O, specCH4, specOH, 3.02000E+05, 2.1, 17422.0, T, X, dWdT, dWdrhok);

  if (REACTION[45])
    add_to_dW_fw_2r2p(specCH3, specO, specCH2O, specH, 6.80000E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[46])
    add_to_dW_fw_2r2p(specCH2O, specH, specHCO, specH2, 5.0E+13, 0.0, 3991.0, T, X, dWdT, dWdrhok);

  if (REACTION[47])
    add_to_dW_fw_2r2p(specCH2O, specOH, specHCO, specH2O, 1.4E+14, 0.0, 1100.0, T, X, dWdT, dWdrhok);

  if (REACTION[48])
    add_to_dW_fw_2r2p(specHCO, specH, specCO, specH2, 4.00000E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[49])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_2r3p(specHCO, specM, specCO, specH, specM, 1.60000E+14, 0.0, 14700.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[50])
    add_to_dW_fw_2r2p(specCH3, specO2, specCH3O, specO, 3.0E+13, 0.0, 25652.0, T, X, dWdT, dWdrhok);

  if (REACTION[51])
    add_to_dW_fw_2r2p(specCH3O, specH, specCH2O, specH2, 2.00000E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[52])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_2r3p(specCH3O, specM, specCH2O, specH, specM, 2.40000E+13, 0.0, 28812.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[53])
    add_to_dW_fw_2r2p(specHO2, specHO2, specH2O2, specO2, 2.0E+14, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[54])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_2r3p(specH2O2, specM, specOH, specOH, specM, 1.30000E+17, 0.0, 45500.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[55])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specOH, specOH, specM, specH2O2, specM, 9.86000E+14, 0.0, -5070.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[56])
    add_to_dW_fw_2r2p(specH2O2, specOH, specH2O, specHO2, 1.00000E+13, 0.0, 1800.0, T, X, dWdT, dWdrhok);

  if (REACTION[57])
    add_to_dW_fw_2r2p(specH2O, specHO2, specH2O2, specOH, 2.86000E+13, 0.0, 32790.0, T, X, dWdT, dWdrhok);

  if (REACTION[58])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specOH, specH, specM, specH2O, specM, 2.20000E+22, -2.0, 0.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[59])  {
    for (specM = 0; specM < ns; specM++) {
      add_to_dW_fw_3r2p(specH, specH, specM, specH2, specM, 1.80000E+18, -1.0, 0.0, T, X, dWdT, dWdrhok);
    }
  }

  if (REACTION[60])
    add_to_dW_fw_2r2p(specCH3, specOH, specCH2, specH2O, 7.6E+06, 2.0, 5000.0, T, X, dWdT, dWdrhok);

  if (REACTION[61])
    add_to_dW_fw_2r2p(specCH2, specO, specCO, specH2, 3.0E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[62])
    add_to_dW_fw_2r2p(specCH2, specOH, specCH, specH2O, 1.13E+07, 2.0, 3000.0, T, X, dWdT, dWdrhok);

  if (REACTION[63])
    add_to_dW_fw_2r2p(specCH, specO, specCO, specH, 5.7E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[64])
    add_to_dW_fw_2r2p(specCH, specOH, specHCO, specH, 3.0E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[65])
    add_to_dW_fw_2r2p(specCH, specO2, specHCO, specO, 3.3E+13, 0.0, 0.0, T, X, dWdT, dWdrhok);

  if (REACTION[66])  
    add_to_dW_fw_2r2p ( specCH, specCO2, specHCO, specCO,  8.40e+13, 0.0, 200.0, T, X, dWdT, dWdrhok );   

}

