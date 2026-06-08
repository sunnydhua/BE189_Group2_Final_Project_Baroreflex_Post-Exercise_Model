function [t, BRE, MAP_error, P_star, MAP] = Baroreceptor_Calculations(alpha_hyd, alpha_train)

global T Rs
in_LV_sa;   % loads repo baseline values

% =========================
% REST BASELINE FROM REPO
% ==========================
HR_rest   = 1/T;        % bpm
Rs_rest   = Rs;         % mmHg/(L/min)
SV_rest   = 0.065;      % L/beat
Tc_rest   = 37.0;       % deg C
nu_rest   = 1.0;        % normalized viscosity

% ==========================================
% END-OF-EXERCISE STATE, MODERATE EXERCISE
% ==========================================
HR_ex       = 150;          % bpm
Tc_ex       = 39.0;         % deg C
contract_ex = 1.40;         % 40% above rest

% ======================================
% BAROREFLEX SETTINGS
% ======================================
G0 = 6.0;

% Effector sensitivities
k_HR       = 0.90;     % bpm per mmHg
k_contract = 0.012;    % contractility/mmHg
k_Rs       = 0.12;     % resistance/mmHg

% Recovery timing
tau_Tc_base = 15.0;    % min
tau_HR_base = 3.0;     % min
tau_Rs_base = 7.0;     % min

% Venous pooling parameters
k_pool  = 0.35;
tau_pool = 5.0;

% ==========================================
% ATHLETE-SPECIFIC RESTING PHYSIOLOGY
% ==========================================

% Lower resting MAP for trained individuals
MAP_rest_baseline = HR_rest*SV_rest*Rs_rest;
P_star = MAP_rest_baseline - 10*(alpha_train - 1);

% Greater exercise vasodilation in trained athletes
Rs_ex = Rs_rest/(1 + alpha_train);

% ===================
% TIME GRID
% ===================
dt_rec = 0.1;
t = (0:dt_rec:30)';
N = length(t);

% ===========================
% BODY TEMPERATURE RECOVERY
% ===========================

% Dehydrated subjects cool more slowly
Tc = Tc_rest + (Tc_ex - Tc_rest) .* ...
    exp(-t./(tau_Tc_base/alpha_hyd));

alpha_temp = 1 - 0.4*(Tc - Tc_rest)/2;
alpha_temp = max(0.55,min(1.0,alpha_temp));

% ===========================
% BAROREFLEX GAIN
% ===========================

G_t = G0 .* alpha_train .* alpha_hyd .* alpha_temp;
Grel = G_t ./ G0;

% ===========================
% PREALLOCATE ARRAYS
% ===========================

HR              = zeros(N,1);
contract_factor = zeros(N,1);
Rs_t            = zeros(N,1);

SV              = zeros(N,1);
Q               = zeros(N,1);
MAP             = zeros(N,1);

MAP_error       = zeros(N,1);

nu              = zeros(N,1);
r_rel           = zeros(N,1);

V_pool          = zeros(N,1);

% ===========================
% INITIAL CONDITIONS
% ===========================

HR(1)              = HR_ex;
contract_factor(1) = contract_ex;
Rs_t(1)            = Rs_ex;

% Initial post-exercise pooling
V_pool(1) = 0.15*(1/alpha_hyd);

% Venous return factor
VR = alpha_hyd*(1 - V_pool(1));

SV(1) = SV_rest * contract_factor(1) * VR;

Q(1) = SV(1)*HR(1);

MAP(1) = Q(1)*Rs_t(1);

MAP_error(1) = P_star - MAP(1);

nu(1) = nu_rest*(1 - 0.02*(Tc(1)-Tc_rest));

r_rel(1) = ((Rs_rest*nu(1))/(Rs_t(1)*nu_rest))^(1/4);

% ============================================
% BAROREFLEX RECOVERY MODEL
% ============================================

for i = 2:N

    % Response speeds
    tau_HR = tau_HR_base/Grel(i-1);
    tau_Rs = tau_Rs_base/Grel(i-1);

    % Pressure error
    err = P_star - MAP(i-1);

    % Reflex commands

    HR_cmd = HR_rest + Grel(i-1)*k_HR*err;
    HR_cmd = min(HR_ex,max(HR_rest,HR_cmd));

    contract_cmd = 1.0 + Grel(i-1)*k_contract*err;
    contract_cmd = min(contract_ex,max(1.0,contract_cmd));

    Rs_cmd = Rs_rest + Grel(i-1)*k_Rs*err;
    Rs_cmd = min(1.25*Rs_rest,max(Rs_ex,Rs_cmd));

    % First-order response

    HR(i) = HR(i-1) + ...
        (HR_cmd - HR(i-1))*dt_rec/tau_HR;

    contract_factor(i) = contract_factor(i-1) + ...
        (contract_cmd - contract_factor(i-1))*dt_rec/tau_HR;

    Rs_t(i) = Rs_t(i-1) + ...
        (Rs_cmd - Rs_t(i-1))*dt_rec/tau_Rs;

    % Blood viscosity

    nu(i) = nu_rest*(1 - 0.02*(Tc(i)-Tc_rest));

    % Venous pooling dynamics

    pool_target = ...
        (0.25*(Rs_rest/Rs_t(i))) * (1/alpha_hyd);

    V_pool(i) = V_pool(i-1) + ...
        (pool_target - V_pool(i-1))*dt_rec/tau_pool;

    V_pool(i) = max(0,min(0.6,V_pool(i)));

    % Venous return

    VR = alpha_hyd*(1 - k_pool*V_pool(i));

    % Hemodynamics

    SV(i) = SV_rest * contract_factor(i) * VR;

    Q(i) = SV(i)*HR(i);

    MAP(i) = Q(i)*Rs_t(i);

    % Error tracking

    MAP_error(i) = P_star - MAP(i);

    % Relative capillary radius

    r_rel(i) = ...
        ((Rs_rest*nu(i))/(Rs_t(i)*nu_rest))^(1/4);

end

% ==================================
% BAROREFLEX EFFECTIVENESS (BRE)
% ==================================

dHR  = gradient(HR,dt_rec);
dMAP = gradient(MAP,dt_rec);

BRE = nan(size(t));

valid = abs(dMAP) > 0.02;

BRE(valid) = -dHR(valid)./dMAP(valid);

BRE = movmean(BRE,5,'omitnan');



    