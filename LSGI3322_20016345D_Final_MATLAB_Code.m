clc
clear all

%==========================================================================
%======================TANG JUSTIN HAYSE CHI WING G.=======================
%======================20016345D  Year 4 LSGI Student======================
%==========================================================================
%================= LSGI3322 First Intermediate Report =====================
%============== Read Both Navigation and Observation Files ================
%==========================================================================

%==========================================================================
%======================Read RINEX navigation File==========================
%==========================================================================
[icondata,iconcmap] = imread('LSGI3322_20016345D_Pop_up.png');
z1=msgbox({'Input a navigation file.' ; 'Press OK.'}, 'Reminder','custom', icondata,iconcmap);
uiwait(z1);
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

No_Epoch_Nav = 1;
while ischar(Read_NavLine)
     Read_NavLine = fgetl(NavData);
     if (Read_NavLine == -1)
         No_Epoch_Nav = No_Epoch_Nav - 1;
     break;
     end

%==========================Broadcase Orbit - 1=============================
    Nav(No_Epoch_Nav).PRN=str2num(Read_NavLine(1:2));
    Nav(No_Epoch_Nav).Year=str2num(Read_NavLine(4:5))+2000;
    Nav(No_Epoch_Nav).Month=str2num(Read_NavLine(7:8));
    Nav(No_Epoch_Nav).Day=str2num(Read_NavLine(10:11));
    Nav(No_Epoch_Nav).Hour=str2num(Read_NavLine(13:14));
    Nav(No_Epoch_Nav).Minute=str2num(Read_NavLine(16:17));
    Nav(No_Epoch_Nav).Second=str2num(Read_NavLine(18:22));
    Nav(No_Epoch_Nav).Date_Numerical=datenum(Nav(No_Epoch_Nav).Year,Nav(No_Epoch_Nav).Month,Nav(No_Epoch_Nav).Day,Nav(No_Epoch_Nav).Hour,Nav(No_Epoch_Nav).Minute,Nav(No_Epoch_Nav).Second);
    No_of_weekday_In_Nav=weekday(Nav(No_Epoch_Nav).Date_Numerical)-1;  %Start from sunday
    Nav(No_Epoch_Nav).Time_in_GPS=No_of_weekday_In_Nav*60*60*24+Nav(No_Epoch_Nav).Hour*60*60+Nav(No_Epoch_Nav).Minute*60+Nav(No_Epoch_Nav).Second;
    Nav(No_Epoch_Nav).SV_Clock_Bias=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).SV_Clock_drift=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).SV_Clock_drift_rate=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 2=============================
    Nav(No_Epoch_Nav).IODE=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).Crs=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).Delta_N=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).M0=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 3=============================
    Nav(No_Epoch_Nav).Cuc=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).e=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).Cus=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).sqrt_a=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 4=============================
    Nav(No_Epoch_Nav).Toe_time=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).Cic=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).Omega_0=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).CIS=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 5=============================
    Nav(No_Epoch_Nav).i0=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).Crc=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).Omega=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).Omega_dot=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 6=============================
    Nav(No_Epoch_Nav).IDOT=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).weekNO=str2num(Read_NavLine(42:60));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 7=============================
    Nav(No_Epoch_Nav).SV_accuracy=str2num(Read_NavLine(4:22));
    Nav(No_Epoch_Nav).SV_health=str2num(Read_NavLine(23:41));
    Nav(No_Epoch_Nav).TGD=str2num(Read_NavLine(42:60));
    Nav(No_Epoch_Nav).IODC=str2num(Read_NavLine(61:79));
    Read_NavLine=fgetl(NavData);
%==========================Broadcase Orbit - 8=============================
    Nav(No_Epoch_Nav).Transmission_timeofMessage=str2num(Read_NavLine(4:22));

    No_Epoch_Nav=No_Epoch_Nav+1;
%==================The END of Reading Navigation Data======================
end

%==========================================================================
%======================Read RINEX Observation File=========================
%==========================================================================
[icondata,iconcmap] = imread('LSGI3322_20016345D_Pop_up.png');
z2=msgbox({'Input an Observation File.' ; 'Press OK.'}, 'Reminder','custom',icondata,iconcmap);
uiwait(z2);
[selObs, Obs_pathname] = uigetfile;
Observation_path = strcat(Obs_pathname,selObs);
ObsData = fopen(Observation_path,'r');
fprintf('A pair of coordinate of GPS Receiver is now computing, Pls wait for a moment.\nThank you for your understanding.\n');

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
IndexNo_A = find(contains(Types_of_observation,'P1'));
IndexNo_B = find(contains(Types_of_observation,'P2'));
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
No_of_weekday_In_Obs=weekday(Obs(Epoch_OBS + h).Date_numerical)-1; %Start from sunday
Obs(Epoch_OBS + h).Time_in_GPS=No_of_weekday_In_Obs*60*60*24+Obs(Epoch_OBS + h).Hour*60*60+Obs(Epoch_OBS + h).Minute*60+Obs(Epoch_OBS + h).Second;
end
%======================Find the C1 Observation Data========================
for d = 1:No_of_PRN
Read_ObsLine = fgetl(ObsData);

Obs(Epoch_OBS + d).C1 = str2num(Read_ObsLine(2+16*(IndexNo_C-1):15+16*(IndexNo_C-1)));

Read_ObsLine = fgetl(ObsData);
end
Epoch_OBS = Epoch_OBS + No_of_PRN;    %Total No of observation
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
Approximate_Coor = [Approx_X;Approx_Y;Approx_Z;RxClockError];
%==========Whether both files with the identical PRN number?===============
while (abs(dx(1,1))+abs(dx(2,1))+abs(dx(3,1))+abs(RxClockDiff)) > 10e-8
for i = 1:Epoch_OBS
q = 1;
for b = 1:No_Epoch_Nav
if isequal(Compare_Nav_PRN(1,b),Compare_Obs_PRN(1,i))
RowOfSamePRN(q,i) = b;
q = q + 1;
end
end
No_RowOfSamePRN = sum(RowOfSamePRN~=0);
%===============The nearest GPS time with 4 Hours Vaildation===============
Difference_Minimum = 4*60*60; %4-hours vaildation

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

%==========================Semi-major Axis (a)=============================
semi_major_axis=(Nav(Sate_Row).sqrt_a)^2;
a = semi_major_axis;
%==================Time from Ephemeris Reference Epoch=====================
half_week = 3.5*60*60*24;       %In terms of GPS Time
t_emission = Obs(i).C1 /c;      %Signal transmission time
t = Obs(i).Time_in_GPS - t_emission - RxClockError;

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
%=======================True anomaly (True_Anomaly)========================
True_Anomaly = 2*atan(sqrt((1 + Nav(Sate_Row).e)/(1 - Nav(Sate_Row).e)) * tan(E/2));
%=================Argument of latitude (Coorected_ArgLat)==================
phi = True_Anomaly + Nav(Sate_Row).Omega;
%================Orbit Radius of the GPS satellite position================
r = a * (1 - Nav(Sate_Row).e * cos(E));
delta_r = Nav(Sate_Row).Crs * sin(2*phi) + Nav(Sate_Row).Crc * cos(2*phi);
Corrected_Radius = r + delta_r;
%===================Corrected GPS orbit plane's Inclination================
inclination_i = Nav(Sate_Row).i0 + Nav(Sate_Row).IDOT * tk;
delta_i = Nav(Sate_Row).CIS* sin(2*phi) + Nav(Sate_Row).Cic * cos(2*phi);
Corrected_Inclination = inclination_i + delta_i;
%====================argument of latitude (Corr0ected_phi)=================
delta_phi = Nav(Sate_Row).Cus * sin(2*phi) + Nav(Sate_Row).Cuc * cos(2*phi);
Corrected_ArgLat = phi + delta_phi;
%====================GPS positions in the orbital plane====================
X0 = Corrected_Radius * cos(Corrected_ArgLat);
Y0 = Corrected_Radius * sin(Corrected_ArgLat);
%=======================longitude of ascending node========================
Corrected_omega = Nav(Sate_Row).Omega_0 + (Nav(Sate_Row).Omega_dot - omega_E) *tk - omega_E * Nav(Sate_Row).Toe_time;
%=========Orbital Coordinate System to Terrestrial Coordinate System=======
%======================Earth-centred Earth-fixed Frame=====================
GPS(i,1) = X0 * cos(Corrected_omega)- Y0 * cos(Corrected_Inclination) * sin(Corrected_omega);
GPS(i,2) = X0 * sin(Corrected_omega)+ Y0 * cos(Corrected_Inclination) * cos(Corrected_omega);
GPS(i,3) = Y0 * sin(Corrected_Inclination);
%==================The END of GPS Satellite Positions======================



%==========================================================================
%======================== LSGI3322 Final Report ===========================
%================ Calculation of GPS Satellite Positions ==================
%==========================================================================

%========================== Error Correction ==============================
%============================Earth Rotation================================
WT = omega_E * t_emission;
Rotation_r(1,1) = cos(WT);
Rotation_r(1,2) = sin(WT);
Rotation_r(1,3) = 0;
Rotation_r(2,1) = -sin(WT);
Rotation_r(2,2) = cos(WT);
Rotation_r(2,3) = 0;
Rotation_r(3,1) = 0;
Rotation_r(3,2) = 0;
Rotation_r(3,3) = 1;

%================ Re-calculation of Orbit Coordinate ======================
GPS_Coordinate = GPS * Rotation_r;
%======================== Satellite Clock Error ===========================
SxClockError(i) = c * (Nav(Sate_Row).SV_Clock_Bias + (Nav(Sate_Row).SV_Clock_drift)* tk + (Nav(Sate_Row).SV_Clock_drift_rate)* tk^2);

%==========================================================================
%====================== GPS Receiver Position =============================
%==========================================================================
%===========Geometric distance between satellite and receiver==============
rho(i) = sqrt((Approximate_Coor(1,1)- GPS(i,1))^2 + (Approximate_Coor(2,1)- GPS(i,2))^2 + (Approximate_Coor(3,1)- GPS(i,3))^2);
%======================= Formation of B Matrix ============================
B(i,1) = (Approximate_Coor(1,1)- GPS(i,1))/rho(i);
B(i,2) = (Approximate_Coor(2,1)- GPS(i,2))/rho(i);
B(i,3) = (Approximate_Coor(3,1)- GPS(i,3))/rho(i);
B(i,4) = 1;
%======================= Formation of Function (f)=========================
f(i,1) = Obs(i).C1 - rho(i) + SxClockError(i);
end
%================Calculation of the Correction (dx)========================
dxprevious = RxClockError;
dx = inv(transpose(B) * B) * transpose(B) * f;

Approximate_Coor = Approximate_Coor + dx;
RxClockError = dx(4,1)/c;
RxClockDiff = RxClockError - dxprevious;
iteration = iteration + 1;
end
%=============== The Coordinate of GPS Recevier Position ==================
GPS_Receiver(1,1) = Approximate_Coor(1);
GPS_Receiver(2,1) = Approximate_Coor(2);
GPS_Receiver(3,1) = Approximate_Coor(3);
GPS_Receiver(4,1) = Approximate_Coor(4);

disp('-------------------------------------------------------------------')
fprintf('1. MATLAB Computation Output: \n The GPS Receiver Position (in ECEF Frame) is X: %f, Y: %f, Z: %f.\n', GPS_Receiver(1,1),GPS_Receiver(2,1),GPS_Receiver(3,1));
disp('-------------------------------------------------------------------')

X =  GPS_Receiver(1,1);
Y =  GPS_Receiver(2,1);
Z =  GPS_Receiver(3,1);

%===================== ECEF Frame to WGS84 Datum ==========================
%===========================Basic Constants================================
a = 6378137;            %Physical Constant: Semi-major Axis
f = 1/298.257223563;    %Physical Constant: Ellipsoid Flattening
b = a*(1-f);            %Physical Constant: Semi-minor Axis
%===========================Auxiliary values===============================
P       = sqrt(X^2 + Y^2);
Theta   = atan(Z*a/P*b);
e       = sqrt(((a^2) - (b^2))/a^2); % First eccentricity
e2      = sqrt((a^2 - b^2)/b^2);     % Second eccentricity
%=====================Initial value of Latitude============================
Latitude = atan2(Z,(P*(1-e^2)));
%=============Iteration loop for estimate Latitude(Latitude)===============
for i = 1:10000
    N = a / sqrt(1-e^2*sin(Latitude).^2);               % Prime vertical
    Altitude = (P/cos(Latitude)) - N;                   % Altitude
    Latitude = atan2(Z,(P*(1-e^2*(N/(N+Altitude)))));   % Latitude
end
%====================== Latitude and Longitude ============================
Longitude = (atan2(Y,X))*180/pi;
Latitude   = Latitude*180/pi;

disp('-------------------------------------------------------------------')
fprintf('2. MATLAB Computation Output: \n The GPS Receiver Position (in WGS84 Datum) is X(Long): %f, Y(Lat): %f, Z(Alt): %f.\n', Longitude, Latitude, Altitude);
disp('-------------------------------------------------------------------')
