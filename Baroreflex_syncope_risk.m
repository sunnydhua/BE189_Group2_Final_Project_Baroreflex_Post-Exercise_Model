% filename: baroreflex_syncope_risk.m
%   WHEN is the baroreflex defective in the 30 min after exercise, and how
%   do HYDRATION and TRAINING combine to set syncope (fainting) risk?
%
% WHAT SYNCOPE RISK MEANS HERE (read this!)
%   Syncope is a TRANSIENT event: it happens during the early MAP dip right
%   after you stop, NOT at the 30-min endpoint. Every conscious person's
%   reflex eventually drags MAP back up, so all cases will be model to recover by 30 min
%   Syncope Risk = how DEEP the early dip goes and how
%   LONG MAP lingers in the danger zones before climbing out. A weaker reflex
%   (untrained + dehydrated) dips lower and spends MORE TIME in the red/yellow
%   zones, which is the high-risk signal (even though it eventually escapes).
%
% LITERATURE-BASED ZONES (absolute MAP, mmHg)
%   Red    (high risk)     : MAP < 65   "syncope imminent"  (~SBP < 90, ESC)
%   Yellow (moderate risk) : 65-75       "presyncope risk"   (>5 min = flag)
%   Safe                   : MAP > 75
%   OH line (orthostatic)  : drop of 20 mmHg from P* (ACC/AHA criterion)
%
%   TRAINING TWEAK (new)
%   Fitter individuals have larger venous compliance, so they POOL MORE at
%   exercise cessation (deeper acute dip) but their stronger reflex rebounds
%   faster. Modeled by making the acute pooling p0 depend on alpha_train.
%
%   P* SET TO 93 mmHg (clinically normal resting MAP from the textbook)
%   Resting Rs is re-anchored so the model's resting MAP equals this normal target
%
%   CASES (2x2): trained/untrained x hydrated/dehydrated
%
%   OUTPUTS
%   Figure:     MAP vs time, all 4 cases + tiered danger zones + OH line
%   Console:    start, nadir, time<65, time 65-75, risk label

clear; clc; close all;
global T TS tauS tauD
global Csa Rs RMi RAo dt
in_LV_sa;   % used for cardiac period T at HR_rest

% 1: REST ANCHORS (clinically normal P*)
HR_rest = 1/T;                                 % bpm (from repo)
SV_rest = 0.065;                               % L/beat
P_star  = 93.0;                                % mmHg, normal resting MAP (textbook)
Rs_rest = P_star/(HR_rest*SV_rest);            % anchor rest MAP to P*
Tc_rest = 37.0;                                % deg C

% 2: END-OF-EXERCISE STATE

HR_ex       = 150;        % bpm
Rs_ex       = Rs_rest/2;  % vasodilated
Tc_ex       = 39.0;       % deg C
contract_ex = 1.40;       % peak contractility

% 3: FIXED REFLEX / PHYSICAL CONSTANTS (shared)

G0        = 6.0;     % baseline reflex gain
k_HR      = 0.90;    % bpm per mmHg error (compensatory tachycardia)
k_Rs      = 0.12;    % mmHg/(L/min) per mmHg error (vasoconstriction)
HR_max    = HR_ex;   % ceiling on reflex tachycardia
p0_base       = 0.55;  % baseline preload fraction at cessation (acute pooling)
c_train_pool  = 0.12;  % TWEAK: fitter -> deeper pooling (lower p0)
tau_pre_base  = 4.0;   % min, venous-return restoration
tau_c         = 3.0;   % min, contractility washout
tau_HR_base   = 3.0;   % min, HR recovery
tau_Rs_base   = 7.0;   % min, vascular tone recovery
tau_Tc_base   = 15.0;  % min, body cooling

% 4: LITERATURE DANGER THRESHOLDS (absolute MAP, mmHg)

MAP_red    = 65.0;        % < this = high risk (syncope imminent)
MAP_yellow = 75.0;        % 65-75 = moderate risk (presyncope)
MAP_OH     = P_star - 20; % orthostatic-hypotension line (drop of 20 from P*)
yellow_min = 5.0;         % min in yellow zone that flags moderate risk

% 5: TIME + TEMPERATURE (shared)

dt_rec = 0.1;  t = (0:dt_rec:30)';  N = length(t);
Tc = Tc_rest + (Tc_ex - Tc_rest).*exp(-t/tau_Tc_base);
alpha_temp = max(0.55, min(1.0, 1 - 0.4*(Tc - Tc_rest)/2));
contract = 1 + (contract_ex - 1).*exp(-t/tau_c);   % passive inotropy washout

% 6: DEFINE THE 4 CASES

cases(1) = struct('name','Trained + Hydrated',  'alpha_train',1.25,'d_hyd',0.00,'col',[0.10 0.60 0.20]);
cases(2) = struct('name','Trained + Dehydrated','alpha_train',1.25,'d_hyd',0.12,'col',[0.90 0.55 0.10]);
cases(3) = struct('name','Untrained + Hydrated','alpha_train',0.75,'d_hyd',0.00,'col',[0.10 0.45 0.85]);
cases(4) = struct('name','Untrained + Dehyd.',  'alpha_train',0.75,'d_hyd',0.12,'col',[0.85 0.15 0.15]);
nC = numel(cases);
MAP_all = zeros(N,nC);
fprintf('\n=== Post-Exercise Syncope Risk (P* = %.0f, red<%.0f, yellow %.0f-%.0f) ===\n',...
       P_star, MAP_red, MAP_red, MAP_yellow);

% 7: RUN EACH CASE

for c = 1:nC
   a_train = cases(c).alpha_train;
   d_hyd   = cases(c).d_hyd;
   phi_vol   = 1 - d_hyd;                        % blood-volume / preload ceiling
   alpha_hyd = 1 - 0.6*d_hyd;                    % hydration effect on gain
   p0        = p0_base - c_train_pool*(a_train-1); % TWEAK: trained pools deeper
   Grel = (a_train * alpha_hyd) .* alpha_temp;   % normalized gain over time
   HR=zeros(N,1); Rs_t=zeros(N,1); preload=zeros(N,1); SV=zeros(N,1); MAP=zeros(N,1);
   HR(1)=HR_ex; Rs_t(1)=Rs_ex; preload(1)=p0*phi_vol;
   SV(1)=SV_rest*contract(1)*preload(1);
   MAP(1)=SV(1)*HR(1)*Rs_t(1);
   for i = 2:N
       err = P_star - MAP(i-1);
       HR_cmd = min(HR_max,      max(HR_rest, HR_rest + Grel(i-1)*k_HR*err));
       Rs_cmd = min(1.3*Rs_rest, max(Rs_ex,   Rs_rest + Grel(i-1)*k_Rs*err));
       g = Grel(i-1);
       HR(i)      = HR(i-1)      + (HR_cmd  - HR(i-1))     *dt_rec/(tau_HR_base /g);
       Rs_t(i)    = Rs_t(i-1)    + (Rs_cmd  - Rs_t(i-1))   *dt_rec/(tau_Rs_base /g);
       preload(i) = preload(i-1) + (phi_vol - preload(i-1))*dt_rec/(tau_pre_base/g);
       SV(i)  = SV_rest*contract(i)*preload(i);
       MAP(i) = SV(i)*HR(i)*Rs_t(i);
   end
   MAP_all(:,c) = MAP;
   % --- Tiered syncope metrics ---
   tRed    = sum(MAP <  MAP_red)*dt_rec;                    % min in red
   tYellow = sum(MAP >= MAP_red & MAP < MAP_yellow)*dt_rec; % min in yellow
   [nadir,iN] = min(MAP);
   tRed_all(c)=tRed; tYel_all(c)=tYellow;
   if     tRed > 0,             risk='HIGH (entered red)';
   elseif tYellow > yellow_min, risk='HIGH (>5 min presyncope)';
   elseif tYellow > 0,          risk='MODERATE (brief presyncope)';
   else,                        risk='LOW';
   end
   fprintf('%-22s start %.1f | nadir %.1f @ %.1f min | <65: %.1f m | 65-75: %.1f m | %s\n',...
       cases(c).name, MAP(1), nadir, t(iN), tRed, tYellow, risk);
end

% 8: MAP trajectories with tiered danger zones
figure(1); clf; hold on
ylo = min(MAP_all(:))-3;  yhi = max(P_star, max(MAP_all(:)))+3;
patch([0 30 30 0],[ylo ylo MAP_red MAP_red],              [1.00 0.80 0.80],'EdgeColor','none'); % red
patch([0 30 30 0],[MAP_red MAP_red MAP_yellow MAP_yellow],[1.00 0.95 0.75],'EdgeColor','none'); % yellow
hh = gobjects(nC,1);
for c=1:nC, hh(c)=plot(t,MAP_all(:,c),'-','Color',cases(c).col,'LineWidth',2.2); end
hR = yline(MAP_red,'r-','LineWidth',1.5);
hY = yline(MAP_yellow,'-','Color',[0.85 0.65 0.0],'LineWidth',1.5);
hP = yline(P_star,'k:','LineWidth',1.2);
xlabel('Time post-exercise (min)'); ylabel('MAP (mmHg)')
title('Post-Exercise MAP: Syncope Risk and Recovery Speed')
xlim([0 30]); ylim([ylo yhi]); grid on; box on
legend([hh; hR; hY; hP], ...
[{cases.name}, {'Red <65 (syncope)'}, {'Yellow 65-75 (presyncope)'}, {sprintf('Resting P*=%.0f',P_star)}], ...
  'Location','southeast')
