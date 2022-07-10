%==========================================================================
%======================Read RINEX Navigation File==========================
%==========================================================================
navigation_path='site0900.01n';
Nav_Data=fopen(navigation_path,'r');
%==========Find out the header section "END OF HEADER"=====================
 Read_Nav_Line=fgetl(Nav_Data);
 while ischar(Read_Nav_Line)
 if contains(Read_Nav_Line, 'END OF HEADER');
 break;
 end

 Read_Nav_Line=fgetl(Nav_Data);
 end
%======================Read Navigation Data================================
Epoch_Nav=1;
while ischar(Read_Nav_Line)
Read_Nav_Line=fgetl(Nav_Data);
if (Read_Nav_Line == -1)
Epoch_Nav=Epoch_Nav-1;
break;
end
%==========================CV/ Epoch/ SV CLK===============================
nav(Epoch_Nav).PRN=str2num(Read_Nav_Line(1:2));
nav(Epoch_Nav).Year=str2num(Read_Nav_Line(4:5));
nav(Epoch_Nav).Month=str2num(Read_Nav_Line(7:8));
nav(Epoch_Nav).Day=str2num(Read_Nav_Line(10:11));
nav(Epoch_Nav).Hour=str2num(Read_Nav_Line(13:14));
nav(Epoch_Nav).Minute=str2num(Read_Nav_Line(16:17));
nav(Epoch_Nav).Second=str2num(Read_Nav_Line(18:22));
nav(Epoch_Nav).Date_Numerical=datenum(nav(Epoch_Nav).Year,nav(Epoch_Nav).Month,nav(Epoch_Nav).Day,nav(Epoch_Nav).Hour,nav(Epoch_Nav).Minute,nav(Epoch_Nav).Second);
No_of_weekday_In_nav=weekday(nav(Epoch_Nav).Date_Numerical)-1;  
nav(Epoch_Nav).Time_in_GPS=No_of_weekday_In_nav*60*60*24+nav(Epoch_Nav).Hour*60*60+nav(Epoch_Nav).Minute*60+nav(Epoch_Nav).Second;
nav(Epoch_Nav).SV_Clock_Bias=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).SV_Clock_drift=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).SV_Clock_drift_rate=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 1=============================     
nav(Epoch_Nav).IODE=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).Crs=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).Delta_N=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).M0=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 2============================= 
nav(Epoch_Nav).Cuc=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).e=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).Cus=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).sqrt_a=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 3=============================     
nav(Epoch_Nav).Toe_time=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).Cic=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).Omega_0=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).CIS=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 4=============================     
nav(Epoch_Nav).i0=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).Crc=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).Omega=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).Omega_dot=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 5=============================     
nav(Epoch_Nav).IDOT=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).weekNO=str2num(Read_Nav_Line(42:60));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 6=============================  
nav(Epoch_Nav).SV_accuracy=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).SV_health=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).TGD=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).IODC=str2num(Read_Nav_Line(61:79));
Read_Nav_Line=fgetl(Nav_Data);
%==========================Broadcase Orbit - 7=============================  
nav(Epoch_Nav).Transmission_timeofMessage=str2num(Read_Nav_Line(4:22));
nav(Epoch_Nav).Fit_interval=str2num(Read_Nav_Line(23:41));
nav(Epoch_Nav).Spare_1=str2num(Read_Nav_Line(42:60));
nav(Epoch_Nav).Spare_2=str2num(Read_Nav_Line(61:79));
    
Epoch_Nav=Epoch_Nav+1;
%==================The END of Reading Navigation Data======================
end
