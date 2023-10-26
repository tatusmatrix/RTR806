format compact
clc

IC_m = 7;
IC_phase_deg = 30;
IC_phase = 30 / 180 * pi;

IC_m_cpx = 7 * exp(j*IC_phase);

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
Iin_m_cpx = IC_m_cpx
% provement by calculation
VC_m_cpx = IC_m_cpx * ZC_cpx
VC_m = abs(VC_m_cpx)
VC_phase_degrees = angle(VC_m_cpx) / pi * 180
VRv_m_cpx = VC_m_cpx
IRv_m_cpx = VRv_m_cpx / ZRv_cpx
Iin_m_cpx = IC_m_cpx + IRv_m_cpx
Iin_m_cpx = IC_m_cpx * (ZC_cpx + ZRv_cpx) / ZRv_cpx

VL_m_cpx = Iin_m_cpx * ZL_cpx
VL_m = abs(VL_m_cpx)
VL_phase_degrees = angle(VL_m_cpx) / pi * 180

VR_m_cpx = Iin_m_cpx * ZR_cpx
VR_m = abs(VR_m_cpx)
VR_phase_degrees = angle(VR_m_cpx) / pi * 180

Vin_m_cpx = Iin_m_cpx * Zin_cpx
Vin_m_cpx = VR_m_cpx + VC_m_cpx + VL_m_cpx
Vin_m = abs(Vin_m_cpx)
Vin_phase_degrees = angle(Vin_m_cpx) / pi * 180

Sin_cpx = 1/2 * Vin_m_cpx * (-Iin_m_cpx)'

% LTspice result analysis
T = 1 / f
delta_L = 20.8e-6
delta_fi = delta_L / T * 360