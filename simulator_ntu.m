tic
format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % 光速(m/s)
lambda = c /(1575.42*10^6) ;

Data = 0116 ;
obs_flag = 1 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 900 ;
initial_time = converge_time + 20 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
interpolation_order = 13 ;                      %內插階數

switch Data
    case{ 1210 }
        addpath DG14_1210
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2015 ;     week_num = 1874 ;         month = 12 ;      day = 10 ;      day_of_year = 344 ;
        obs_filename = '_msvgpsuv.txt' ;                                      % GPS的數據檔案
        SAT_filename = 'sat_data_V1A1.csv' ;                               % SAT資料
        MOT_filename = 'motion_V1.csv' ;                               % MOT資料
        
end
        
%-----------讀取接收機資訊------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

[ time_ms , SAT_PRN , SAT_PX , SAT_PY , SAT_PZ , SAT_VX , SAT_VY , SAT_VZ , AZ , EL , RANGE , prL1 , pr_rate , I_delay , T_delay , signal_dB , Ant_AZ , Ant_EL , Range_rate , D_s , D_range , SAT_AX , SAT_AY , SAT_AZ ] = ...
    textread( SAT_filename , '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 2 ) ;%headerlines:跳過前兩行

%[ time_ms , SAT_PRN , SAT_PX , SAT_PY , SAT_PZ , SAT_VX , SAT_VY , SAT_VZ , AZ , EL , RANGE , prL1 , pr_rate , I_delay , T_delay , signal_dB , Ant_AZ , Ant_EL , Range_rate , D_s , D_range , SAT_AX , SAT_AY , SAT_AZ ] = ...
%    textread( SAT_filename , '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 2 ) ;%headerlines:跳過前兩行


%{
time(1:300) = [] ;
prn(1:300) = [] ;
SV_x(1:300) = [] ;
SV_y(1:300) = [] ;
SV_z(1:300) = [] ;
ADR(1:300) = [] ;
pr(1:300) = [] ;
pr_rate(1:300) = [] ;
SV_x_dot(1:300) = [] ;
SV_y_dot(1:300) = [] ;
SV_z_dot(1:300) = [] ;
    %}
%-----------initialize--------------------------
SVMax = 32 ;
runtime = round( max( time ) ) - round ( min( time ) ) +1;
%runtime = 1200 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
nsat_Mtrix = EDOP ;
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
wl_count = 0 ;     %計算GPS/BDS權重疊代次數
chg = 1 ; % record the change of satellite amounts
PRr = zeros( 32,2 ) ;
Carrier = zeros( 32,2 ) ;
Smooth_Count = zeros( 32,1 ) ;
Smooth_Code = zeros( 32,1 ) ;

%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8


%----------------猜測初始位置--------------------------
llh_0 = [ 25.039535 121.517567 3 ] ;       %景福門
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
rec_bias = 0 ;
first_postime = 0;
bar = waitbar(0,'Please wait...') ;

for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    rec_bias_satpos = 0 ;
    id = [] ;
    Rho0 = [] ;
    Rho0_pr = [] ;
    CNR = [] ;
    EL = [] ;
    AZI = [] ;
    PR = zeros( 32,1 ) ;
    PRr( :,1 ) = PRr( :,2 ) ;
    PRr( :,2 ) = 0 ;
    Carrier( :,1 ) = Carrier( :,2 ) ;
    Carrier( :,2 ) = 0 ;
    SVxyz = zeros( 1,3 ) ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    
    while judge == 1
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
            if prn( count ) ~= 193  && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999 %&&prn( count ) ~= 6&&prn( count ) ~= 24
                
                SVxyz(1,1) = SV_x(count) ;
                SVxyz(1,2) = SV_y(count) ;
                SVxyz(1,3) = SV_z(count) ;
                
                PR( prn(count) , 1 ) = pr( count ) ;
                PRr( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier( prn(count) , 2 ) = ADR( count ) ;          %* lambda ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    id( j , 1) = prn( count ) ;
                    [ EL( j,1) , AZI( j,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , SVxyz ) ;
                    Rho0_pr( j , 1 ) = pr( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && Smooth_Count( prn(count),1 ) > converge_time && Carrier(prn(count),2) ~= 0 && Carrier(prn(count),1) ~= 0
                    id( j , 1) = prn( count ) ;
                    [ EL( j,1) , AZI( j,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , SVxyz ) ;
                    j = j+1 ;
                end
                
            end
            count = count + 1 ;
            
        else
            count_times = time( count ) ;
            judge = 0 ;
        end
        
        if count>length( time )
            count = count - 1 ;
            break
        end
        
    end
    
    %----------------每次(每一秒)定位衛星的顆數---------------
    nsat = length( id ) ;
    nsat_Mtrix( i ) = nsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %-------------------計算covariance----------------------------
        CNR = zeros(nsat,1) ;
        [ sat_var_elevation , sat_var_CNR , sat_var_SIGMA , sat_var_CandE ] = weightingfunc( EL , CNR ) ;
        %-------------------選擇權重矩陣------------------------------
        W = diag( ( 1./( sat_var_elevation ) ).^2 ) ;
        %W = eye( nsat ) ;
        
        %----------------選擇觀測量------------------
        if obs_flag == 1
            Rho0 = Rho0_pr ;
        elseif obs_flag == 2
            [ Rho0 , Smooth_Code , Smooth_Count ] = Doppler_Smoothed_Code( PRr , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        elseif obs_flag == 3
            [ Rho0 , Smooth_Code , Smooth_Count ] = Carrier_Smoothed_Code( Carrier , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        end
        
        Rr = Rho0 ;
        Delta_Mtrix(1:4) = 1 ;
        
        while norm( Delta_Mtrix(1:3) ) > 0.01
            
            llh_0 = xyz2llh( rec_pos_0 ) ;
            lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
            lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
            height_rec = llh_0( 3 ) ;                  % unit : meter
            
            Tropo = zeros( nsat , 1 ) ;
            Iono = zeros( nsat , 1 ) ;
            
            %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------
            % 北斗的SCBb包含sat_bias,relative與group_delay
            Rhoc = Rho0 - Iono - Tropo + SCBg + relative_IGS ;
            
            %---------------------------修正地球自轉--------------------------------
            [ Xs , Rr ] = Fix_Earth_Rotation( sat , rec_pos_0 ) ;
            
            %-----------------------最小平方法定位-------------------------
            if nsat > 3     %   單系統至少需要4顆衛星來定位
                [ Delta_Mtrix , Delta_Rho ] = LeastSquare( Xs , rec_pos_0 , Rhoc , W , nsat , rec_bias ) ;
                
                rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                rec_bias = rec_bias + Delta_Mtrix(4) ;
                rec_pos_0 = rec_pos ;
            else
                break
            end
            
        end
        
        if nsat < 4
            Xxyz( i , : ) = rec_pos_0 ;
        else
            Xxyz( i , : ) = rec_pos ;
            if first_postime ==0
                first_postime = i
            end
        end
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
        rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置
        
    end
    
    %----------------時間+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    if i > 1 && nsat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    
    nsat_prior = nsat ; %   紀錄上一秒衛星顆數
    
end

close( bar ) ;
switch obs_flag
    case { 0 }
        obs = 'EPR' ;
    case { 1 }
        obs = 'PR' ;
    case{ 2 }
        obs = 'DSC' ;
    case{ 3 }
        obs = 'CSC' ;
end
%---------------------資料取樣時間----------------------
sampling_time_end = runtime ;
effTimeArray = first_postime : sampling_time_end ;
effTimeArray_error = round(length(effTimeArray)/100) -1 + first_postime : sampling_time_end ;

%---------------------顯示參數設定---------------------
[ Data_start_end_obs ] = { Data sampling_time_end obs  }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray_error,:),1 ) ;      %計算STD
%rec_pos_act = lla2ecef( [25.149442 , 121.777873 , 73.6] , 'WGS84' );
Xenu = zeros( runtime , 3 ) ;
MSE_enu = zeros(1,3);

for i = 1 : runtime
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
 
end
Xenu = Xenu - solid_tide_error ;
for i = first_postime : runtime
   MSE_enu = MSE_enu + Xenu( i , : ).^2;
end
%--------------------計算STD & RMS------------------------
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;
std = sqrt(length(effTimeArray) \ MSE_enu) ;
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
STD = [ std_enu' ; std_en' ]


%rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
%rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ;
%RMS = [ rms_enu' ; rms_en' ] ;

%error = [ std' , rms' ]
%Error = [ STD , RMS ]

figure( 1 ) ;
plot( effTimeArray , PDOP( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS)' )
print( '-dpng',  'PDOP_GPS' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 2 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS)' )
print( '-dpng',  'Number of visible satellites_GPS' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 3 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
%axis( [-4 , 4 ,-4 , 4 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-20,20] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI


toc