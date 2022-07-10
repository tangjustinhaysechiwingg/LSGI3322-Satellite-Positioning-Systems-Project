clc
clear all
%==========================================================================
%================= LSGI3322 First Intermediate Report =====================
%============== Read Both Navigation and Observation Files ================
%==========================================================================

%==========================================================================
%======================Read RINEX navigation File==========================
%==========================================================================
f1=msgbox({'Input a navigation file.' ; 'Press OK.'}, 'Reminder','warn');
uiwait(f1);
[selNav, Nav_pathname] = uigetfile;
Navigation_path = strcat(Nav_pathname,selNav);
NavData = fopen(Navigation_path,'r');
%==========Find out the header section "END OF HEADER"=====================
Read_NavLine = fgetl(NavData);
while ischar(Read_NavLine)
if contains(Read_NavLine, 'END OF HEADER')
break;
end
Read_NavLine = fgetl(NavData);
end
%======================Read RINEX navigation File==========================
Epoch_NAV = 1;
while ischar(Read_NavLine)
Read_NavLine = fgetl(NavData);
if (Read_NavLine == -1)
Epoch_NAV = Epoch_NAV - 1;
break;
end
%==========================Broadcase Orbit - 1=============================    
Nav(Epoch_NAV).PRN=str2num(Read_NavLine(1:2));
Nav(Epoch_NAV).Year=str2num(Read_NavLine(4:5))+2000;
Nav(Epoch_NAV).Month=str2num(Read_NavLine(7:8));
Nav(Epoch_NAV).Day=str2num(Read_NavLine(10:11));
Nav(Epoch_NAV).Hour=str2num(Read_NavLine(13:14));
Nav(Epoch_NAV).Minute=str2num(Read_NavLine(16:17));
Nav(Epoch_NAV).Second=str2num(Read_NavLine(18:22));
Nav(Epoch_NAV).Date_Numerical=datenum(Nav(Epoch_NAV).Year,Nav(Epoch_NAV).Month,Nav(Epoch_NAV).Day,Nav(Epoch_NAV).Hour,Nav(Epoch_NAV).Minute,Nav(Epoch_NAV).Second);
No_of_weekday_In_Nav=weekday(Nav(Epoch_NAV).Date_Numerical)-1;  
Nav(Epoch_NAV).Time_in_GPS=No_of_weekday_In_Nav*60*60*24+Nav(Epoch_NAV).Hour*60*60+Nav(Epoch_NAV).Minute*60+Nav(Epoch_NAV).Second;
Nav(Epoch_NAV).SV_Clock_Bias=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).SV_Clock_drift=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).SV_Clock_drift_rate=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 2=============================  
Nav(Epoch_NAV).IODE=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).Crs=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).Delta_N=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).M0=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 3=============================  
Nav(Epoch_NAV).Cuc=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).e=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).Cus=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).sqrt_a=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 4=============================  
Nav(Epoch_NAV).Toe_time=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).Cic=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).Omega_0=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).CIS=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 5=============================  
Nav(Epoch_NAV).i0=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).Crc=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).Omega=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).Omega_dot=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 6=============================  
Nav(Epoch_NAV).IDOT=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).weekNO=str2num(Read_NavLine(42:60));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 7=============================  
Nav(Epoch_NAV).SV_accuracy=str2num(Read_NavLine(4:22));
Nav(Epoch_NAV).SV_health=str2num(Read_NavLine(23:41));
Nav(Epoch_NAV).TGD=str2num(Read_NavLine(42:60));
Nav(Epoch_NAV).IODC=str2num(Read_NavLine(61:79));
Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 8=============================  
Nav(Epoch_NAV).Transmission_timeofMessage=str2num(Read_NavLine(4:22));
    
Epoch_NAV=Epoch_NAV+1;
%==================The END of Reading Navigation Data======================
end
%==========================================================================
%======================Read RINEX Observation File=========================
%==========================================================================
f2=msgbox({'Input an observation file.' ; 'Press OK.'}, 'Reminder','warn');
uiwait(f2);
[selObs, Obs_pathname] = uigetfile;
Observation_path = strcat(Obs_pathname,selObs);
ObsData = fopen(Observation_path,'r');
fprintf('Please Find the Satellite Position on "Workplace".');
%============Find the Wording: "'APPROX POSITION XYZ'"=====================
Read_Obs_Approx_XYZ_Line = fgetl(ObsData);
while ischar(Read_Obs_Approx_XYZ_Line)
if contains(Read_Obs_Approx_XYZ_Line, 'APPROX POSITION XYZ')
break;
end
Read_Obs_Approx_XYZ_Line = fgetl(ObsData);
end
Approx_X = str2num(Read_Obs_Approx_XYZ_Line(2:14));
Approx_Y = str2num(Read_Obs_Approx_XYZ_Line(16:28));
Approx_Z = str2num(Read_Obs_Approx_XYZ_Line(30:42));
%============Find the Wording: "'# / TYPES OF OBSERV'"=====================
Read_Obs_C_Line = fgetl(ObsData);
while ischar(Read_Obs_C_Line)
if contains(Read_Obs_C_Line, '# / TYPES OF OBSERV')
break;
end
Read_Obs_C_Line = fgetl(ObsData);
end
No_of_TypesOfObservation = str2num(Read_Obs_C_Line(4:6));
%=====Find the C1 Position in the '# / TYPES OF OBSERV'====================
for p = 1:No_of_TypesOfObservation
    Types_of_observation{p} = Read_Obs_C_Line(4+6*p:6+6*p);   
end
IndexNo_C = find(contains(Types_of_observation,'C1'));
%============Find the Wording: "'TIME OF FIRST OBS'"=======================
Read_ObsLine = fgetl(ObsData);
while ischar(Read_ObsLine)
if contains(Read_ObsLine, 'TIME OF FIRST OBS')
Obs_Time_of_FirstObs = string(Read_ObsLine(5:6));
break;
end
Read_ObsLine = fgetl(ObsData);
end
%==================Find the Wording: "END OF HEADER"=======================
Read_ObsLine = fgetl(ObsData);
while ischar(Read_ObsLine)
if contains(Read_ObsLine, 'END OF HEADER')
break;
end
Read_ObsLine = fgetl(ObsData);
end
%======================Read Observation Data===============================  
Epoch_OBS = 0;
while ischar(Read_ObsLine)
Read_ObsLine = fgetl(ObsData);
if (Read_ObsLine == -1)
break;
end
%=====================Find the Observation Time============================
if isequal(string(Read_ObsLine(2:3)),Obs_Time_of_FirstObs)
No_of_PRN = str2num(Read_ObsLine(31:32));
%==================Find the total Number of Satellites=====================
for h = 1: No_of_PRN
Obs(Epoch_OBS + h).PRN=str2num(Read_ObsLine(31+3*h:32+3*h));
Obs(Epoch_OBS + h).Year=str2num(Read_ObsLine(1:3))+2000;
Obs(Epoch_OBS + h).Month=str2num(Read_ObsLine(5:6));
Obs(Epoch_OBS + h).Day=str2num(Read_ObsLine(7:9));
Obs(Epoch_OBS + h).Hour=str2num(Read_ObsLine(11:12));
Obs(Epoch_OBS + h).Minute=str2num(Read_ObsLine(14:15));
Obs(Epoch_OBS + h).Second=str2num(Read_ObsLine(17:26));
Obs(Epoch_OBS + h).Epoch_Flag=str2num(Read_ObsLine(28:29));
Obs(Epoch_OBS + h).Date_numerical=datenum(Obs(Epoch_OBS + h).Year,Obs(Epoch_OBS + h).Month,Obs(Epoch_OBS + h).Day,Obs(Epoch_OBS + h).Hour,Obs(Epoch_OBS + h).Minute,Obs(Epoch_OBS + h).Second);
No_of_weekday_In_Obs=weekday(Obs(Epoch_OBS + h).Date_numerical)-1;
Obs(Epoch_OBS + h).Time_in_GPS=No_of_weekday_In_Obs*60*60*24+Obs(Epoch_OBS + h).Hour*60*60+Obs(Epoch_OBS + h).Minute*60+Obs(Epoch_OBS + h).Second;
end
%======================Find the C1 Observation Data========================
for d = 1:No_of_PRN
Read_ObsLine = fgetl(ObsData);
Obs(Epoch_OBS + d).C1 = str2num(Read_ObsLine(2+16*(IndexNo_C-1):15+16*(IndexNo_C-1)));  
Read_ObsLine = fgetl(ObsData);
end
Epoch_OBS = Epoch_OBS + No_of_PRN;
end      
%==================The END of Reading Observation Data=====================
end
%========Extract PRN Dada through Navigation and Observation file==========
Compare_Nav_PRN=[Nav(:).PRN];
Compare_Obs_PRN=[Obs(:).PRN];
%============GPS Time of Navigation and Observation file===================
Compare_Nav_GPSTime=[Nav(:).Time_in_GPS];
Compare_Obs_GPSTime=[Obs(:).Time_in_GPS];
%==========Match the GPS Time of Navigation and Observation file===========
RxClockError = 0.01; 
RxClockDiff = 1;
iteration=0;
dx = [0.01;0.01;0.01;RxClockError];
Approx_Coordinate = [Approx_X;Approx_Y;Approx_Z;RxClockError];
%==========Whether both files with the identical PRN number?===============
for i = 1:Epoch_OBS
q = 1; 
for b = 1:Epoch_NAV
if isequal(Compare_Nav_PRN(1,b),Compare_Obs_PRN(1,i))
RowOfSamePRN(q,i) = b; 
q = q + 1; 
end  
end
No_RowOfSamePRN = sum(RowOfSamePRN~=0);
%===============The nearest GPS time with 4 Hours Vaildation===============
Difference_Minimum = 4*60*60;  
for m = 1:No_RowOfSamePRN(i)
Difference_GPSTime = abs(Compare_Obs_GPSTime(1,i) - Compare_Nav_GPSTime(1,RowOfSamePRN(m,i)));
if (Difference_GPSTime < Difference_Minimum)
Difference_Minimum = Difference_GPSTime;
Sate_Row = RowOfSamePRN(m,i);
end            
end
%==========================================================================
%================= LSGI3322 Second Intermediate Report ====================
%================ Calculation of GPS Satellite Positions ==================
%==========================================================================
  
%===========================Basic Constants================================
c=299792458;                %Physcial Constant: Speed of Light (m/s)
g=9.80665;                  %Physcial Constant: Accelera
G=6.67259e-11;              %Physcial Constant: Constant of Gravity (Nm^2/kg^2)
GM=3.986005e+14;            %Spaceflight Constant: GM(Earth) (m^3/s^2)
omega_E=7.2921151467e-05;   %Earth rotation rate (rad/s)
%==========================Keplerian Elements==============================
%==========================Semi-major Axis (a)=============================
semi_major_axis=(Nav(Sate_Row).sqrt_a)^2;
a = semi_major_axis;
%==================Time from Ephemeris Reference Epoch=====================
half_week = 3.5*60*60*24;       %In terms of GPS Time
Transmission_Time = Obs(i).C1 /c;      %Signal transmission time
t = Obs(i).Time_in_GPS - Transmission_Time - RxClockError;

tk = (t - Nav(Sate_Row).Toe_time);
if (tk > half_week)
tk = tk - (2*half_week);
elseif (tk < -half_week)
tk = tk + (2*half_week);
end
%==========================Mean motion (n)=================================
n0 = sqrt (GM / (a^3));
n = n0 + Nav(Sate_Row).Delta_N;
%==========================Mean anomaly (M)================================
M = Nav(Sate_Row).M0 + (n * tk);
%==========================Eccentric anomaly (E)===========================
n = 1;
E = 1;
E0 = M;
while n <= 30 && abs(E0 - E) > 10^(-14)
E_new = M +(Nav(Sate_Row).e) * sin(E);
E = E0;
E0 = E_new;
n = n + 1;
end
%=========================True anomaly (True_Anomaly)=============================
True_Anomaly = 2*atan(sqrt((1 + Nav(Sate_Row).e)/(1 - Nav(Sate_Row).e)) * tan(E/2));
%==================Argument of latitude (Coorected_ArgLat)====================
phi = True_Anomaly + Nav(Sate_Row).Omega;
%================Orbit Radius of the GPS satellite position================
r = a * (1 - Nav(Sate_Row).e * cos(E));
delta_r = Nav(Sate_Row).Crs * sin(2*phi) + Nav(Sate_Row).Crc * cos(2*phi);
Orbit_Radius = r + delta_r;
%===================Corrected GPS orbit plane's Inclination================
inclination_i = Nav(Sate_Row).i0 + Nav(Sate_Row).IDOT * tk;
delta_i = Nav(Sate_Row).CIS* sin(2*phi) + Nav(Sate_Row).Cic * cos(2*phi);
Corrected_Inclination = inclination_i + delta_i;
%====================argument of latitude (Corr0ected_phi)=================
delta_phi = Nav(Sate_Row).Cus * sin(2*phi) + Nav(Sate_Row).Cuc * cos(2*phi);
Coorected_ArgLat = phi + delta_phi;
%====================GPS positions in the orbital plane====================
X0 = Orbit_Radius * cos(Coorected_ArgLat);
Y0 = Orbit_Radius * sin(Coorected_ArgLat);
%=======================longitude of ascending node========================
Omega_K = Nav(Sate_Row).Omega_0 + (Nav(Sate_Row).Omega_dot - omega_E) *tk - omega_E * Nav(Sate_Row).Toe_time;
%=========Orbital Coordinate System to Terrestrial Coordinate System=======
%======================Earth-centred Earth-fixed Frame=====================
Satellite_Position(i,1) = X0 * cos(Omega_K)- Y0 * cos(Corrected_Inclination) * sin(Omega_K);
Satellite_Position(i,2) = X0 * sin(Omega_K)+ Y0 * cos(Corrected_Inclination) * cos(Omega_K);
Satellite_Position(i,3) = Y0 * sin(Corrected_Inclination);
%==================The END of GPS Satellite Positions======================
end
