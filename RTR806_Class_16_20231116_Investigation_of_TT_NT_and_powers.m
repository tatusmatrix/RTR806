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

% Thevenin's and Norton's theorems
% Thevenin's
V3m_cpx = 0;
%ZRload_cpx_TT = 1.e9;

YL1R1_cpx = 1 / (ZR1_cpx + ZL1_cpx)
YC1_cpx = 1 / ZC1_cpx
YC2_cpx = 1 / ZC2_cpx
%YL3Rload_cpx_TT = 1 / (ZL3_cpx + ZRload_cpx_TT)
YL3Rload_cpx_TT = 0;

Y_cpx = [(YL1R1_cpx + YC1_cpx + YC2_cpx) -(YC2_cpx)
                              -(YC2_cpx)  (YC2_cpx + YL3Rload_cpx_TT)]

Im_cpx = [ Vs1m_cpx * YL1R1_cpx  ;  Is2m_cpx ]

Vnodesm_cpx = Y_cpx \ Im_cpx

V2m_cpx = Vnodesm_cpx(2);

Vmm_cpx = V2m_cpx
Vnm_cpx = V3m_cpx
Vmnocm_cpx = Vmm_cpx - Vnm_cpx
Vstm_cpx = Vmnocm_cpx

Vstm = abs(Vstm_cpx)
fi_of_Vst_deg = angle(Vstm_cpx) / pi * 180

% Norton's
ZRload_cpx_NT = 0;
Z_cpx = [ (ZR1_cpx + ZC1_cpx + ZL1_cpx) -(ZC1_cpx)
                             -(ZC1_cpx)  (ZC1_cpx+ZC2_cpx+ZL3_cpx+ZRload_cpx_NT)]
Vm_cpx = [ Vs1m_cpx ; -Is2m_cpx*(ZL3_cpx+ZRload_cpx_NT)]

Imeshm_cpx = Z_cpx \ Vm_cpx

Imnscm_cpx = Imeshm_cpx(2) + Is2m_cpx
Isnm_cpx = Imnscm_cpx

Isnm = abs(Isnm_cpx)
fi_of_Isn_deg = angle(Isnm_cpx) / pi * 180


% Zt_cpx and Zn_cpx
Zeq1_cpx = ZR1_cpx + ZL1_cpx
Zeq2_cpx = Zeq1_cpx * ZC1_cpx / (Zeq1_cpx + ZC1_cpx)
Zeq_cpx = Zeq2_cpx + ZC2_cpx + ZL3_cpx
Zt = Zeq_cpx
Zt = Vstm_cpx / Isnm_cpx

Rt = real(Zt)
Xt = imag(Zt)
if Xt > 0
  Lt = Xt / w
  Cload_eq = 1 / (Xt * w)
else
  Ct = abs(1 / (Xt * w))
  Lload_eq = abs(Xt / w)
end

% simulation of powers
%Rload_ = Rload
Rload_ = Rt*0.001 : 0.01 : Rt*10;
Cload_ = Cload_eq
Xload_ = -1 / (w * Cload_)
Zload_cpx_ = Rload_ + j*Xload_;

Rsource = Rt
Lsource = Lt
Xsource = w * Lsource
Zsource_cpx = Rsource + j*Xsource

Iload_ttm_cpx = Vstm_cpx./(Zsource_cpx + Zload_cpx_);
VRload_ttm_cpx = Rload_.*Iload_ttm_cpx;
SRload_tt_cpx = 1/2 * VRload_ttm_cpx.*conj(Iload_ttm_cpx);
PRload_tt_cpx = real(SRload_tt_cpx);

plot(Rload_, PRload_tt_cpx), grid on, hold on

Cload_ = 0.999*Cload_eq
Xload_ = -1 / (w * Cload_)
Zload_cpx_ = Rload_ + j*Xload_;

Iload_ttm_cpx = Vstm_cpx./(Zsource_cpx + Zload_cpx_);
VRload_ttm_cpx = Rload_.*Iload_ttm_cpx;
SRload_tt_cpx = 1/2 * VRload_ttm_cpx.*conj(Iload_ttm_cpx);
PRload_tt_cpx = real(SRload_tt_cpx);

plot(Rload_, PRload_tt_cpx), hold on

Cload_ = 0.998*Cload_eq
Xload_ = -1 / (w * Cload_)
Zload_cpx_ = Rload_ + j*Xload_;

Iload_ttm_cpx = Vstm_cpx./(Zsource_cpx + Zload_cpx_);
VRload_ttm_cpx = Rload_.*Iload_ttm_cpx;
SRload_tt_cpx = 1/2 * VRload_ttm_cpx.*conj(Iload_ttm_cpx);
PRload_tt_cpx = real(SRload_tt_cpx);

plot(Rload_, PRload_tt_cpx), hold on
