% from Lab 2 Controlled CV
%filename:  LV_sa.m
clear all %clear all variables
clf       %and figures
global T TS tauS tauD;
global Csa Rs RMi RAo dt CHECK PLA;
in_LV_sa  %initialize
for klok=1:klokmax
  t=klok*dt;
  PLV_old=PLV;
  Psa_old=Psa;
  CLV_old=CLV;
  CLV=CV_now(t,CLVS,CLVD);
  %find self-consistent 
  %valve states and pressures:
  set_SMi_SAo
  %store in arrays for future plotting:
  t_plot(klok)=t;
  CLV_plot(klok)=CLV;
  PLV_plot(klok)=PLV;
  Psa_plot(klok)=Psa;
  VLV_plot(klok)=CLV*PLV+VLVd;
  Vsa_plot(klok)=Csa*Psa+Vsad;
  QMi_plot(klok)=SMi*(PLA-PLV)/RMi;
  QAo_plot(klok)=SAo*(PLV-Psa)/RAo;
  Qs_plot(klok)=Psa/Rs;
  SMi_plot(klok)=SMi;
  SAo_plot(klok)=SAo;
end

%plot results:
figure(1)
subplot(3,1,1), plot(t_plot,CLV_plot)
  ylabel('C_{LV} (L/mmHg)'); title('Task 2: Csa = 0.001156 L/mmHg')
subplot(3,1,2), plot(t_plot,PLV_plot,t_plot,Psa_plot)
  ylabel('Pressure (mmHg)'); legend('P_{LV}','P_{sa}')
subplot(3,1,3), plot(t_plot,QMi_plot,t_plot,QAo_plot,t_plot,Qs_plot)
  ylabel('Flow (L/min)'); xlabel('t (min)'); legend('Q_{Mi}','Q_{Ao}','Q_s')

%left ventricular pressure-volume loop
figure(2)
plot(VLV_plot,PLV_plot)
  xlabel('V_{LV} (L)'); ylabel('P_{LV} (mmHg)'); title('LV Pressure-Volume Loop')

%systemic arterial pressure-volume ``loop''
figure(3)
plot(Vsa_plot,Psa_plot)
  xlabel('V_{sa} (L)'); ylabel('P_{sa} (mmHg)'); title('Systemic Arterial P-V')

save('LV_sa.mat')
