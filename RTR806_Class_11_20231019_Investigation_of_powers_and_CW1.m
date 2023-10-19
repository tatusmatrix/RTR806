format compact
clc

Ic_m = 7;
Ic_phase_deg = 30;
Ic_phase = 30 / 180 * pi;

Ic_m_cpx = 7 * exp(j*Ic_phase);

f = 10.e3;

R = 100;     Rv = 100.e6;
C = 5.e-6;
L = 5.e-3;

w = 2*pi*f
ZR_cpx = R
ZRv_cpx = Rv
XC = -1/(w*C)
XL = w*L
ZC_cpx = j * XC
ZC_cpx = -j / (w*C)
ZL_cpx = j * XL

Zparallel_cpx = ZC_cpx * ZRv_cpx / (ZC_cpx + ZRv_cpx)
Yparallel_cpx = 1/ZC_cpx + 1/ZRv_cpx
Zparallel_cpx = 1/Yparallel_cpx

Zin_cpx = ZR_cpx + Zparallel_cpx + ZL_cpx
Lin = imag(Zin_cpx) / w

% logical reason
Iin_m_cpx = Ic_m_cpx

Vin_m_cpx = Iin_m_cpx * Zin_cpx

Sin_cpx = 1/2 * Vin_m_cpx * (-Iin_m_cpx)'