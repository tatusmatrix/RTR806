clear
clc
format compact

f = 1.e3;
w = 2*pi*f;

Vs1m = 15;
fi_of_Vs1_deg = -30;
fi_of_Vs1_rad = fi_of_Vs1_deg / 180 * pi;
Vs1m_cpx = Vs1m * exp(j*fi_of_Vs1_rad)

Is2m = 3;
fi_of_Is2_deg = 45;
fi_of_Is2_rad = fi_of_Is2_deg / 180 * pi;
Is2m_cpx = Is2m * exp(j*fi_of_Is2_rad)

R1 = 3;
C1 = 20.e-6; C2 = 50.e-6;
L1 = 10.e-3; L3 = 20.e-3;
Rload = 50;

ZR1_cpx = R1;
XC1 = -1/(w*C1);
ZC1_cpx = j*XC1;
XC2 = -1/(w*C2);
ZC2_cpx = j*XC2;
XL1 = w*L1;
ZL1_cpx = j*XL1;
XL3 = w*L3;
ZL3_cpx = j*XL3;
ZRload_cpx = Rload;

% Mesh method:

Z_cpx = [ (ZR1_cpx + ZC1_cpx + ZL1_cpx) -(ZC1_cpx)
                             -(ZC1_cpx)  (ZC1_cpx+ZC2_cpx+ZL3_cpx+ZRload_cpx)]
Vm_cpx = [ Vs1m_cpx ; -Is2m_cpx*(ZL3_cpx+ZRload_cpx)]

Imeshm_cpx = Z_cpx \ Vm_cpx

IRloadm_cpx = Imeshm_cpx(2) + Is2m_cpx
VRloadm_cpx = ZRload_cpx * IRloadm_cpx

IRloadm = abs(IRloadm_cpx)
fi_of_IRload_deg = angle(IRloadm_cpx) / pi * 180

VRloadm = abs(VRloadm_cpx)
fi_of_VRload_deg = angle(VRloadm_cpx) / pi * 180

% Nodal analysis
V3m_cpx = 0;

YL1R1_cpx = 1 / (ZR1_cpx + ZL1_cpx)
YC1_cpx = 1 / ZC1_cpx
YC2_cpx = 1 / ZC2_cpx
YL3Rload_cpx = 1 / (ZL3_cpx + ZRload_cpx)

Y_cpx = [(YL1R1_cpx + YC1_cpx + YC2_cpx) -(YC2_cpx)
                              -(YC2_cpx)  (YC2_cpx + YL3Rload_cpx)]

Im_cpx = [ Vs1m_cpx * YL1R1_cpx  ;  Is2m_cpx ]

Vnodesm_cpx = Y_cpx \ Im_cpx

V2m_cpx = Vnodesm_cpx(2);

IRloadm_cpx = (V2m_cpx - V3m_cpx) / (ZL3_cpx + ZRload_cpx)
VRloadm_cpx = ZRload_cpx * IRloadm_cpx
