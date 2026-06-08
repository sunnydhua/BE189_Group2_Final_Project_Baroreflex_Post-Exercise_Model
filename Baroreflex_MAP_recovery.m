% filename: baroreflex_MAP_recovery.m  
%	After the full 30-min post-exercise recovery, is the individual's MAP
%	HYPOTENSIVE, NORMAL, or HYPERTENSIVE
%	Healthy reflex: returns to normal
%	Chronotropic incompetence: HR can't compensate and linger hypotensive
%	Sympathetic overactivity: set-point reset up and hypertensive
%	Autonomic failure: no vasoconstriction and stay hypotensive
%
% BP BANDS (AHA readings converted to MAP via MAP = (2*DBP + SBP)/3)
%     120/80 -> 93.3 mmHg  (top of normal)
%     130/80 -> 96.7 mmHg  (stage-1 hypertension)
%      90/60  -> 70.0 mmHg  (hypotension)
%   => Hypotensive < 70 | Normal 70-93.3 | Elevated 93.3-96.7 | Hypertensive >= 96.7
%
%   NOTE ON SCALING
%   Resting Rs is anchored so resting MAP = MAP_rest_healthy (a clinically
%   NORMAL value), making endpoints directly comparable to the AHA bands.
%
% OUTPUTS
%   Figure:     MAP trajectories per case with shaded BP bands
%   Console:    final MAP + label per case

clear; clc; close all;
global T TS tauS tauD
global Csa Rs RMi RAo dt
in_LV_sa;   % repo baseline (used for cardiac period T -> HR_rest)

% 1: REST ANCHORS (clinically normal baseline)

HR_rest = 1/T;             % bpm (from repo)
SV_rest = 0.065;           % L/beat
MAP_rest_healthy = 88;     % mmHg, mid-normal resting MAP target
Rs_rest = MAP_rest_healthy/(HR_rest*SV_rest);  % anchor rest MAP to normal

% 2: END-OF-EXERCISE STATE

HR_ex = 150; Rs_ex = Rs_rest/2; contract_ex = 1.40; HR_max = 150;
Tc_rest = 37.0; Tc_ex = 39.0;

% 3: FIXED REFLEX / PHYSICAL CONSTANTS

k_HR = 0.90; k_Rs = 0.12; p0 = 0.55;
tau_HR0 = 3.0; tau_Rs0 = 7.0; tau_pre0 = 4.0; tau_c = 3.0; tau_Tc = 15.0;

% 4: BP CLASSIFICATION BANDS (MAP, mmHg)

MAP_hypo = 70.0;    % below = hypotensive
MAP_norm = 93.3;    % top of normal (120/80)
MAP_htn  = 96.7;    % hypertensive (130/80)

% 5: TIME + TEMPERATURE

dt_rec = 0.1; t = (0:dt_rec:30)'; N = numel(t);
Tc = Tc_rest + (Tc_ex - Tc_rest).*exp(-t/tau_Tc);
alpha_temp = max(0.55, min(1.0, 1 - 0.4*(Tc - Tc_rest)/2));
contract = 1 + (contract_ex - 1).*exp(-t/tau_c);   % passive inotropy washout

% 6) CASES (each is a baroreflex phenotype; unique pattern of blood pressure reflex response)
%   delta_set : baroreflex set-point shift (+ = resetting upward)
%   gain      : overall reflex gain scale (autonomic strength)
%   chrono    : chronotropic competence (1 normal, low = HR can't rise)
%   vaso      : vasoconstrictor competence (1 normal, low = can't restore tone)
%   phi       : recoverable blood-volume / preload ceiling (1 = full)
cases(1)=struct('name','Healthy reflex',         'delta',0.00,'gain',1.20,'chrono',1.00,'vaso',1.00,'phi',1.00);
cases(2)=struct('name','Chronotropic incompet.', 'delta',0.00,'gain',0.80,'chrono',0.35,'vaso',0.95,'phi',0.85);
cases(3)=struct('name','Sympathetic overactivity','delta',0.18,'gain',1.00,'chrono',1.00,'vaso',1.00,'phi',1.00);
cases(4)=struct('name','Autonomic failure',      'delta',0.00,'gain',0.40,'chrono',0.60,'vaso',0.60,'phi',0.90);
nC = numel(cases);
MAP_all = zeros(N,nC);

% 7: RUN EACH CASE

fprintf('\n=== Post-30-min BP Classification (bands: hypo<%.0f | normal | htn>=%.1f) ===\n',...
       MAP_hypo, MAP_htn);
for c = 1:nC
   P_set   = MAP_rest_healthy*(1 + cases(c).delta);          % defended set-point
   Grel    = cases(c).gain .* alpha_temp;                    % effective gain
   HR_ceil = HR_rest + cases(c).chrono*(HR_max - HR_rest);   % reachable HR
   Rs_ceil = cases(c).vaso*1.3*Rs_rest;                      % reachable tone
   phi     = cases(c).phi;
   HR=zeros(N,1); Rs_t=zeros(N,1); pre=zeros(N,1); MAP=zeros(N,1);
   HR(1)=HR_ex; Rs_t(1)=Rs_ex; pre(1)=p0*phi;
   MAP(1)=SV_rest*contract(1)*pre(1)*HR(1)*Rs_t(1);
   for i = 2:N
       err = P_set - MAP(i-1);
       HR_cmd = min(HR_ceil, max(HR_rest, HR_rest + Grel(i-1)*k_HR*err));
       Rs_cmd = min(Rs_ceil, max(Rs_ex,  Rs_rest + Grel(i-1)*k_Rs*err));
       g = Grel(i-1);
       HR(i)   = HR(i-1)   + (HR_cmd - HR(i-1)) *dt_rec/(tau_HR0 /g);
       Rs_t(i) = Rs_t(i-1) + (Rs_cmd - Rs_t(i-1))*dt_rec/(tau_Rs0 /g);
       pre(i)  = pre(i-1)  + (phi    - pre(i-1)) *dt_rec/(tau_pre0/g);
       MAP(i)  = SV_rest*contract(i)*pre(i)*HR(i)*Rs_t(i);
   end
   MAP_all(:,c) = MAP;
   % Classify the 30-min endpoint
   mf = MAP(end);
   if     mf <  MAP_hypo, label='HYPOTENSIVE';
   elseif mf <  MAP_norm, label='NORMAL';
   elseif mf <  MAP_htn,  label='ELEVATED';
   else,                  label='HYPERTENSIVE';
   end
   fprintf('%-26s final MAP = %5.1f mmHg  -> %s\n', cases(c).name, mf, label);
end

% 8: Trajectories with BP bands
figure(1); clf; hold on
xlims=[0 30]; ylo=min(MAP_all(:))-4; yhi=max(MAP_all(:))+4;
patch([0 30 30 0],[ylo ylo MAP_hypo MAP_hypo],[1 0.85 0.85],'EdgeColor','none'); % hypo
patch([0 30 30 0],[MAP_htn MAP_htn yhi yhi],  [0.85 0.88 1],'EdgeColor','none'); % htn
cols=[0.10 0.60 0.20; 0.90 0.55 0.10; 0.10 0.45 0.85; 0.85 0.15 0.15];
hh=gobjects(nC,1);
for c=1:nC, hh(c)=plot(t,MAP_all(:,c),'-','Color',cols(c,:),'LineWidth',2.2); end
yline(MAP_norm,'k:','LineWidth',1.2); yline(MAP_htn,'b--','LineWidth',1.2);
yline(MAP_hypo,'r--','LineWidth',1.2);
xlabel('Time post-exercise (min)'); ylabel('MAP (mmHg)')
title('30-min Recovery: Hypotensive vs. Normal vs. Hypertensive Endpoints')
xlim(xlims); ylim([ylo yhi]); grid on; box on
legend(hh,{cases.name},'Location','east')
