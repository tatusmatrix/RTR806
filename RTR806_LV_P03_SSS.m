format compact, clear, clc

% uzdotais lielums - us(t)=Us*cos(w*t+fi_us) - harmoniskā laika funkcija

Us = 10                                             % svārstību amplitūda, V
w = 200*pi                                          % leņķiskā frekvence, rad/s
fi_us_deg = 0                                       % _sākuma_ fāze (laika momentā t = 0s), °
fi_us_rad = fi_us_deg / 180 * pi                    % -"-, radiāni

R = 10;                                             % Ω
L = 1/(10*pi)                                       % H
C = 1/(1000*pi)                                     % F
Ra = 10;                                            % Ω
La = 1/(20*pi)                                      % H

% NB! j vērtību nemainām, jo pēc neklusējuma tā ir imaginārā vienība Matlab'ā
Us_cpx = Us * exp(j*fi_us_rad)                      % izmantota eksponenciālā forma, V
Us_cpx = Us * (cos(fi_us_rad) + j*sin(fi_us_rad))   % -"- trigonometriskā -"-, V
Us = abs(Us_cpx)                                    % modulis - svārstību amplitūda, V
Us_real = real(Us_cpx)                              % reālā daļa, V
Us_imag = imag(Us_cpx)                              % iamginārā daļa, V
Us = sqrt(Us_real^2 + Us_imag^2)                    % modulis - svārstību amplitūda, V
fi_us_rad = angle(Us_cpx)                           % _sākuma_ fāze (laika momentā t = 0s), radiāni
fi_us_deg = fi_us_rad / pi * 180                    % -"-, °

% Z_cpx - kompleksā impedance, Ω
% Z_cpx  = R + j * X
% R - aktīvā pretestība, Ω
% X - reaktīvā pretestība, Ω
ZR_cpx = R + j * 0                                  % Ω; ZR_cpx = R - sastāv tikai no aktīvās pretestības
ZL_cpx = 0 + j * w * L                              % Ω; ZL_cpx = j * w * L - sastāv tikai no XL = w * L, XL > 0
XL = w * L                                          % Ω
ZL_cpx = 0 + j * XL
ZC_cpx = 0 - j / (w * C)                            % Ω; ZC_cpx = -j / (w * C) - sastāv tikai no XC = -1 / (w * C), XC < 0
XC = -1 / (w * C)                                   % Ω
ZC_cpx = 0 + j * XC

ZRa_cpx = Ra + j * 0
ZLa_cpx = 0 + j * w * La
XLa = w * La
ZLa_cpx = 0 + j * XLa

Zeq_cpx = ZR_cpx + ZL_cpx + ZC_cpx                  % virknes slēgums - drīkst saskaitīt Z_cpx
Zeqa_cpx = ZRa_cpx + ZLa_cpx                        % -"-
Zin_cpx = Zeqa_cpx * Zeq_cpx / (Zeqa_cpx + Zeq_cpx) % paralēlais slēgums (zari satur tikai pasīvus elementus)
Yeq_cpx = 1 / Zeq_cpx                               % S vai 1/Ω
Yeqa_cpx = 1 / Zeqa_cpx                             % S vai 1/Ω
Yin_cpx = Yeqa_cpx + Yeq_cpx                        % paralēlais slēgums - drīkst saskaitīt Y_cpx (zari satur tikai pasīvus elementus)
Zin_cpx = 1 / Yin_cpx

I_cpx = Us_cpx / Zeq_cpx
Ia_cpx = Us_cpx / Zeqa_cpx
Is_cpx = - I_cpx - Ia_cpx                            % I_cpx noteikšana, pielietojot KStL A mezglam
Is_cpx = -Us_cpx / Zin_cpx                           % I_cpx noteikšana, izmantojot strāvu caur Zin_cpx, tāpēc prētējā zīme

UR_cpx = I_cpx * ZR_cpx
UL_cpx = I_cpx * ZL_cpx
UC_cpx = I_cpx * ZC_cpx
URa_cpx = Ia_cpx * ZRa_cpx
ULa_cpx = Ia_cpx * ZLa_cpx

% S_cpx - kompleksā jauda, [|S_cpx|] = VA
% S_cpx = P + j * Q
% P - aktīvā jauda, W
% Q - reaktīvā jauda, VAr
% NB! ir jāņem strāvas kompleksi saistītais lielums (fāze ar prētējo zīmi => imaginārā daļa ar prētējo zīmi)
% reizinātājs 1 / 2, ja aprēķinā tiek izmantotas strāvu un spriegumu amplitūdu vērtības
SR_cpx = 1/2 * UR_cpx * I_cpx'                      % SR_cpx = PR - sastāv tikai no aktīvās jaudas
delta_fi_R_rad = angle(UR_cpx) - angle(I_cpx)       % fāzes nobīde starp spriegumu un strāvu, rad
delta_fi_R_deg = delta_fi_R_rad / pi * 180          % -"-, °; R gadījumā - 0°
PR = 1/2 * abs(UR_cpx) * abs(I_cpx) * cos(delta_fi_R_rad)
PR = 1/2 * abs(UR_cpx)^2 / R
PR = 1/2 * abs(I_cpx)^2 * R

SL_cpx = 1/2 * UL_cpx * I_cpx'                      % SL_cpx = j * QL - sastāv tikai no reaktīvās jaudas, QL > 0
delta_fi_L_rad = angle(UL_cpx) - angle(I_cpx)
delta_fi_L_deg = delta_fi_L_rad / pi * 180          % -"-, °; L gadījumā - 90°
QL = 1/2 * abs(UL_cpx) * abs(I_cpx) * sin(delta_fi_L_rad)
QL = 1/2 * abs(UL_cpx)^2 / XL
QL = 1/2 * abs(I_cpx)^2 * XL

SC_cpx = 1/2 * UC_cpx * I_cpx'                      % SC_cpx = j * QC - sastāv tikai no reaktīvās jaudas, QC < 0
delta_fi_C_rad = angle(UC_cpx) - angle(I_cpx)
delta_fi_C_deg = delta_fi_C_rad / pi * 180          % -"-, °; C gadījumā - -90°
QC = 1/2 * abs(UC_cpx) * abs(I_cpx) * sin(delta_fi_C_rad)
QC = 1/2 * abs(UC_cpx)^2 / XC
QC = 1/2 * abs(I_cpx)^2 * XC

SRa_cpx = 1/2 * URa_cpx * Ia_cpx'
delta_fi_Ra_rad = angle(URa_cpx) - angle(Ia_cpx)
delta_fi_Ra_deg = delta_fi_Ra_rad / pi * 180
PRa = 1/2 * abs(URa_cpx) * abs(Ia_cpx) * cos(delta_fi_Ra_rad)
PRa = 1/2 * abs(URa_cpx)^2 / Ra
PRa = 1/2 * abs(Ia_cpx)^2 * Ra

SLa_cpx = 1/2 * ULa_cpx * Ia_cpx'
delta_fi_La_rad = angle(ULa_cpx) - angle(Ia_cpx)
delta_fi_La_deg = delta_fi_La_rad / pi * 180
QLa = 1/2 * abs(ULa_cpx) * abs(Ia_cpx) * sin(delta_fi_La_rad)
QLa = 1/2 * abs(ULa_cpx)^2 / XLa
QLa = 1/2 * abs(Ia_cpx)^2 * XLa

Ss_cpx = 1/2 * Us_cpx * Is_cpx'

% Teledžena teorēma
TT_cpx = SR_cpx + SL_cpx + SC_cpx + SRa_cpx + SLa_cpx + Ss_cpx
TT_cpx = Ss_cpx + (PR + PRa) + j * (QL + QC + QLa)