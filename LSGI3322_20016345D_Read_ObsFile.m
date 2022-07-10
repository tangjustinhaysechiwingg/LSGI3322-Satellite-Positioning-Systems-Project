%==========================================================================
%======================Read RINEX Observation File=========================
%==========================================================================
Observation_path='site0900.01o';
ObsData=fopen(Observation_path,'r');

%======================Find the number of Satellites=======================
Read_Obs_C_Line=fgetl(ObsData);
while ischar(Read_Obs_C_Line)
if contains(Read_Obs_C_Line, '# / TYPES OF OBSERV');
break;
end

Read_Obs_C_Line=fgetl(ObsData);
end

No_of_Satellite=str2num(Read_Obs_C_Line(4:6));
%===============Find the C Position Index of Satellites====================
for i=1:No_of_Satellite
Types_of_observation{i}=Read_Obs_C_Line(4+6*i:6+6*i);   
end

C_position_index = find(contains(Types_of_observation,'C1'));

%===============Find the Wording: "TIME OF FIRST OBS"======================
Read_ObsLine=fgetl(ObsData);
while ischar(Read_ObsLine)
if contains(Read_ObsLine, 'TIME OF FIRST OBS');
Obs_Time_of_first_obs=string(Read_ObsLine(5:6));
break;
end

Read_ObsLine=fgetl(ObsData);
end
%===============Find the Wording: "END OF HEADER"==========================
Read_ObsLine=fgetl(ObsData);
while ischar(Read_ObsLine)
if contains(Read_ObsLine, 'END OF HEADER');
break;
end

Read_ObsLine=fgetl(ObsData);
end
%======================Read Navigation Data================================  
Epoth_Obs=0;
while ischar(Read_ObsLine)
Read_ObsLine=fgetl(ObsData);
if (Read_ObsLine == -1)
Epoth_Obs=Epoth_Obs-1;
break;
end

if isequal(string(Read_ObsLine(2:3)),Obs_Time_of_first_obs)
   No_of_PRN=str2num(Read_ObsLine(31:32));

for j=1: No_of_PRN 
Obs(Epoth_Obs + j).PRN=str2num(Read_ObsLine(31+3*j:32+3*j));
Obs(Epoth_Obs + j).Year=str2num(Read_ObsLine(1:3))+2000;
Obs(Epoth_Obs + j).Month=str2num(Read_ObsLine(5:6));
Obs(Epoth_Obs + j).Day=str2num(Read_ObsLine(7:9));
Obs(Epoth_Obs + j).Hour=str2num(Read_ObsLine(11:12));
Obs(Epoth_Obs + j).Minute=str2num(Read_ObsLine(14:15));
Obs(Epoth_Obs + j).Second=str2num(Read_ObsLine(17:26));
Obs(Epoth_Obs + j).Epoch_Flag=str2num(Read_ObsLine(28:29));
Obs(Epoth_Obs + j).Date_numerical=datenum(Obs(Epoth_Obs + j).Year,Obs(Epoth_Obs + j).Month,Obs(Epoth_Obs + j).Day,Obs(Epoth_Obs + j).Hour,Obs(Epoth_Obs + j).Minute,Obs(Epoth_Obs + j).Second);
No_of_weekday_In_Obs=weekday(Obs(Epoth_Obs + j).Date_numerical)-1;
Obs(Epoth_Obs + j).Time_in_GPS=No_of_weekday_In_Obs*60*60*24+Obs(Epoth_Obs + j).Hour*60*60+Obs(Epoth_Obs + j).Minute*60+Obs(Epoth_Obs + j).Second;
end
     
for k=1:No_of_PRN
Read_ObsLine=fgetl(ObsData);
Obs(Epoth_Obs + k).C1 = str2num(Read_ObsLine(2+16*C_position_index:14+16*C_position_index));  
Read_ObsLine=fgetl(ObsData);
end

Epoth_Obs=Epoth_Obs + No_of_PRN;
end 
%==================The END of Reading Observation Data=====================
end
