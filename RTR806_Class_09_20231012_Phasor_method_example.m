format compact
clc

Iin_abs = 5;
Iin_phase = 45;
Iin_cpx = Iin_abs*exp(j*Iin_phase/180*pi)

w = 1.e3;
f = w/2/pi

R1 = 3;
R2 = 15;
R3 = 8;
L1 = 10.e-3; XL1 =     w*L1
L3 = 16.e-3; XL3 =     w*L3
C2 = 25.e-6; XC2 = -1/(w*C2)

ZR1_cpx = R1 + j*0
ZR2_cpx = R2 + j*0
ZR3_cpx = R3 + j*0
ZL1_cpx = 0 + j*XL1
ZL3_cpx = 0 + j*XL3
ZC2_cpx = 0 + j*XC2

Zeq1_cpx = ZR1_cpx + ZL1_cpx
Zeq1_abs = abs(Zeq1_cpx)
Zeq1_angle = angle(Zeq1_cpx)/pi*180

Zeq2_cpx = ZC2_cpx * ZR2_cpx / (ZC2_cpx + ZR2_cpx)
Zeq2_abs = abs(Zeq2_cpx)
Zeq2_angle = angle(Zeq2_cpx)/pi*180

Zeq3_cpx = ZR3_cpx * ZL3_cpx / (ZR3_cpx + ZL3_cpx)
Zeq3_abs = abs(Zeq3_cpx)
Zeq3_angle = angle(Zeq3_cpx)/pi*180

Zin_cpx = Zeq1_cpx + Zeq2_cpx + Zeq3_cpx
Zin_abs = abs(Zin_cpx)
Zin_angle = angle(Zin_cpx)/pi*180


Vin_cpx = Iin_cpx * Zin_cpx
Vin_abs = abs(Vin_cpx)
Vin_phase = angle(Vin_cpx)/pi*180

VL1_cpx = Iin_cpx * ZL1_cpx
VR1_cpx = Iin_cpx * ZR1_cpx
VC2_cpx = Iin_cpx * Zeq2_cpx
VR2_cpx = VC2_cpx
VR3_cpx = Iin_cpx * Zeq3_cpx
VL3_cpx = VR3_cpx

IL1_cpx = Iin_cpx
IR1_cpx = Iin_cpx
IC2_cpx = VC2_cpx / ZC2_cpx
IR2_cpx = VR2_cpx / ZR2_cpx
IR3_cpx = VR3_cpx / ZR3_cpx
IL3_cpx = VL3_cpx / ZL3_cpx

SL1_cpx = 1/2 * VL1_cpx * IL1_cpx'
SR1_cpx = 1/2 * VR1_cpx * IR1_cpx'
SC2_cpx = 1/2 * VC2_cpx * IC2_cpx'
SR2_cpx = 1/2 * VR2_cpx * IR2_cpx'
SR3_cpx = 1/2 * VR3_cpx * IR3_cpx'
SL3_cpx = 1/2 * VL3_cpx * IL3_cpx'

Sin_cpx = 1/2 * Vin_cpx * (-Iin_cpx)'

Tellegen_s_Theorem = SL1_cpx + SR1_cpx + SC2_cpx + SR2_cpx + SR3_cpx + SL3_cpx + Sin_cpx