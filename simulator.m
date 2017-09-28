tic
format long g
clc ;
clear ;
close all ;
addpath Functions


Data = 1207  ;
no_solid_tide = Data ;
obs_flag = 1 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 900 ;
initial_time = converge_time + 20 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
G_B_condition = 1.05 ;
interpolation_order = 13 ;                      %內插階數
Remark = '' ;

switch Data
    case { 1207 }
        addpath Data_1207 ;
        %rec_pos_act = lla2ecef( [ 25.019672 121.540897 58.5 ] , 'WGS84' ) ;
        year = 2013 ;      month = 12 ;     day = 6 ;       day_of_year = 341 ;
        sampling_time_start = 20 ;                       %sampling_time_end = 10991 ;                  solid_filename = 'solid_1117.txt' ;
        obs_filename = 'COM3_2016-12-07_14.57.43.13o' ; 
                             nav_filename = 'COM3_2016-12-07_14.57.43.13n' ;
end

c = 299792458 ;                                 % 光速(m/s)
L1_freq = 1575.42*10^6 ;                % L1 band frequency(Hz)
lambda_G = c /(1575.42*10^6) ;

%-----------------------讀取接收機資料--------------------------------------------------
[XYZ_station,obs,observablesHeader,measurementsInterval]=readRinexBS(obs_filename);

week_num = obs( : , 1 ) ;
time = obs( : , 2) ;
prn = obs(: , 4) ;
pr = obs(: , 5);
ADR = obs( : , 6) ;
pr_rate = obs(: , 7) ;
cnr0 = obs(: , 8) ;
%-----------------------Preparation-----------------------------------------------
runtime = round( max( time ) ) - round ( min( time ) ) + 1 ;
runtime = 3600 ;
%-----------initialize--------------------------
GMax = 32 ;
BMax = 15 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
EDOP_G=EDOP;                        NDOP_G=EDOP;     VDOP_G=EDOP;     HDOP_G=EDOP;     PDOP_G=EDOP;     GDOP_G=EDOP;
EDOP_B=EDOP;                         NDOP_B=EDOP;     VDOP_B=EDOP;     HDOP_B=EDOP;      PDOP_B=EDOP;     GDOP_B =EDOP;
nGsat_Mtrix = EDOP ;              nBsat_Mtrix = EDOP ;                              nsat_Mtrix = EDOP ;
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
sum0 = 0 ;
chg = 1 ;                                               % record the change of satellite amounts
Rrg = zeros( 32,2 ) ;                                  Rrb = zeros( 15,2 ) ;
Carrier_G = zeros( 32,2 ) ;                       Carrier_B = zeros( 15,2 ) ;
count_G = zeros( 32,1 ) ;                         count_B = zeros( 15,1 ) ;
Delta_CP_G = zeros( runtime , 32 ) ;      Delta_CP_B = zeros( runtime , 15 ) ;
S_G = zeros( 32,1 ) ;                                  S_B = zeros( 15,1 ) ;
SV_Center = zeros( runtime,3 ) ;
ADR_G = zeros( runtime,32 ) ;               ADR_B = zeros( runtime,15 ) ;
Delta_set = zeros( runtime,5 ) ;
Sat2s = zeros( 46,2 ) ;
Cx_set = zeros( runtime,3 ) ;
gpst = zeros( runtime,1 ) ;
RecBias = zeros( runtime, 2 ) ;     % column 1 : GPS     , column 2 : BDS
Delta_total = cell( runtime,4 ) ;
Cx = zeros( 1,3 ) ;
dx = zeros( 1,3 ) ;
GL_count = 1 ;      G_count = 1 ;       L_count = 1 ;
rec_bias_satpos = 0 ;
GP_id_set = zeros( runtime,32 ) ;
BD_id_set = zeros( runtime,15 ) ;
%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week+10800 ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;


%--------------------------------讀取廣播星歷資料----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage

%---------------solid earth tide error----------------------------------
if Data ~= no_solid_tide
    [ time_solid , n_solid , e_solid , u_solid ] = textread( solid_filename ,'%f %f %f %f','delimiter',' ','headerlines',2); %headerlines:跳過前兩行
end
%----------------猜測初始位置--------------------------
%llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
%llh_0 = [ 22.285340 114.161677 3 ] ;       %香港摩天輪
%llh_0 = [-6.367849  106.833466 3 ] ;        %cibg清真寺
%llh_0 = [35, 115 3 ] ;
%rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
rec_pos_0 =  [-2144855.42, 4397605.31, 4078049.85];
%rec_pos_0 = [-1837002.6304  6065627.3599  -716183.2716];
rec_bias_G = 0 ;
rec_bias = 0 ;

Data

runtime_GEL=zeros(32,runtime);
run_trueRange=zeros(32,runtime);
run_ADR=zeros(32,runtime);
run_recbias = zeros(runtime) ;
first_P = 0 ;       %   用來記錄第一次成功定位的時刻數
bar = waitbar(0,'Please wait...');



for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    
    hel_count=0;
    wl_count = 0 ;                                      %計算GPS/BDS權重疊代次數
    GP_id = [] ;
    Rho0_G = [] ;
    Rho0_Gpr = [] ;
    CNR_G = [] ;
    EL_G = [] ;
    AZ_G = [] ;
    Rg = zeros( 32,1 ) ;
    Rrg( :,1 ) = Rrg( :,2 ) ;
    Rrg( :,2 ) = 0 ;
    Carrier_G( :,1 ) = Carrier_G( :,2 ) ;
    Carrier_G( :,2 ) = 0 ;
    SVxyz = zeros( 1,3 ) ;
    AEL_G = [] ;            AAZ_G = [] ;
    ARr_G = [] ;            ARho0_G = [] ;
    ASCBg = [] ;            Asat_G = [] ;
    AGP_id = [] ;            ACNR_G = [] ;
    Atime_tx = [] ;
    
    
    stEL_G = [] ;
    m = 1 ;
    l = 1 ;
    s = 1 ;
    first_in = 0 ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    k = 1 ;                           %每一秒都把各顆BDS的資料依序抓進來
    
    while judge == 1
        
        
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
            if prn( count ) ~= 1193%&&prn( count ) ~= 1013&&prn( count ) ~= 1023 %&& all( AX( prn(count)-1000,: ) ) == 1 && AS( prn(count)-1000,1 ) ~= 999999.999999...
                    prn( count ) = prn( count ) -1000 ;
                
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_G( prn(count) , 2 ) = ADR( count )*lambda_G;
                
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    GP_id( j , 1) = prn( count ) ;
                    Rho0_Gpr( j , 1 ) = pr( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                    GP_id( j , 1) = prn( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
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
    nGsat = length( GP_id ) ;
    nsat = nGsat ;
    
    
    %-----------------------------------------------------------------------------------------------
    if nsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %----------------選擇觀測量------------------
        if obs_flag == 1
            Rho0_G = Rho0_Gpr ;
        elseif obs_flag == 2
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
        elseif obs_flag == 3
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
        end
        
        Rr_G = Rho0_G ;
        
        
        Tg = zeros( nGsat , 1 ) ;
        ig = zeros( nGsat , 1 ) ;
        
        %--------------------得到transmission time ( time_tx )----------------------
        [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
            satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
        %------------使用 time_tx 計算精密星歷位置-----------------
        %--------------initialize----------------
        sat_G = XS ;
        SCBg = c .* dtS ;
        
        %-------------------計算仰角----------------------------------
        for a= 1 : nGsat
            [ EL_G( a , 1 ) , AZ_G( a,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , sat_G( a , : ) ) ;
            %if EL_G( a , 1 ) > smaller_than_elevation
                AEL_G( l , 1) = EL_G( a , 1 ) ;
                AAZ_G( l , 1) = AZ_G( a , 1 ) ;
                ARr_G( l , 1) = Rr_G( a , 1 ) ;
                ARho0_G( l , 1) = ARr_G( l , 1) ;
                ASCBg( l , 1) = SCBg( a , 1 ) ;
                Asat_G( l , :) = sat_G(a , :) ;
                AGP_id( l , 1) = GP_id(a , 1) ;
                Atime_tx( l , 1) = time_tx( a , 1) ;
                runtime_GEL(AGP_id(l,1),i)=AEL_G( l , 1);
               ACNR_G( l , 1) = CNR_G( a , 1);
                l = l + 1 ;
            %else
             %   stEL_G(s,1)= GP_id( a , 1 ) ;
            %    s = s + 1 ;
           % end
        end
        
        AnGsat = length( AGP_id ) ;
        Ansat = AnGsat  ;
        %-------------------計算covariance----------------------------
        EL =  AEL_G  ;
        CNR = ACNR_G ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        
        %-------------------選擇權重矩陣------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( Ansat ) ;
        if AnGsat ~=0
            
            Delta_Mtrix( 1:4 ) = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 0.01
                
                llh_0 = xyz2llh( rec_pos_0 ) ;
                lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
                lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
                height_rec = llh_0( 3 ) ;                  % unit : meter
                
                Tg = zeros( AnGsat , 1 ) ;
                ig = zeros( AnGsat , 1 ) ;
                
                %--------------------得到transmission time ( time_tx )----------------------
                [ XS , dtS , XS_tx , VS_tx , Atime_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                    satellite_positions( r_gpst , ARr_G , AGP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
                
                sat_G = XS ;
                SCBg = c .* dtS ;

                %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------
                % 北斗的SCBb包含sat_bias,relative與group_delay
                Rhoc_G = ARho0_G ;
                Rhoc = [ Rhoc_G ] ;
                
                
                
                %---------------------------修正地球自轉--------------------------------
                [ Xs_G , ARr_G ] = Fix_Earth_Rotation( sat_G , rec_pos_0 ) ;
                Xs = [ Xs_G ] ;
                
                %-----------------------最小平方法定位-------------------------
                if Ansat > 3
                    
                    [ Delta_Mtrix , Delta_Rho ] = LeastSquare( Xs , rec_pos_0 , Rhoc , W , Ansat , rec_bias ) ;
                    
                    rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias = rec_bias + Delta_Mtrix(4) ;
                    rec_pos_0 = rec_pos ;
                    
                    
                    
                else
                    break
                end
                if hel_count >= 100
                    i
                    norm( Delta_Mtrix(1:3) )
                    AGP_id
                    ABD_id
                    
                    break
                end
            end
            
            if Ansat < 4
                Xxyz( i , : ) = rec_pos_0 ;
            else
                Xxyz( i , : ) = rec_pos ;
            end
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置
            %--------------將低仰角衛星之CSC資料清除----------
            %{
            if obs_flag > 1
                for q = 1 : size(stEL_G)
                    S_G(stEL_G(q),1) = 0 ;
                    count_G(stEL_G(q),1) = 0 ;
                end
                for q = 1 : size(stEL_B)
                    S_G(stEL_B(q),1) = 0 ;
                    count_B(stEL_B(q),1) = 0 ;
                end
            end
            %}
        end
    end
    %------------紀錄carrier phase&true range-----------------
    for r=1:AnGsat
        run_ADR(GP_id(r),i) = Carrier_G( AGP_id(r) , 2 ) + ig(r , 1) - Tg(r , 1)- rec_bias ;
        
        %run_trueRange(GP_id(r),i) = sqrt(sat_G(r,1)^2+sat_G(r,2)^2+sat_G(r,3)^2) ;
        run_trueRange(GP_id(r),i) = ARr_G(r) ;
    end
    run_recbias(i) = rec_bias;
    %----------------時間+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    %---------------------計算固體潮汐誤差(solid tide)----------------
    if Data ~= no_solid_tide
        e_solid_error = interp1( time_solid , e_solid , time_day_UTC ) ;
        n_solid_error = interp1( time_solid , n_solid , time_day_UTC ) ;
        u_solid_error = interp1( time_solid , u_solid , time_day_UTC ) ;
        solid_tide_error_temp = [ e_solid_error , n_solid_error , u_solid_error ] ;
        solid_tide_error( i , : ) = solid_tide_error_temp ;
    end
    
    if i > 1 && Ansat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    
    nsat_prior = nsat ; %   紀錄上一秒衛星顆數
    nGsat_Mtrix( i ) = AnGsat ;
    nsat_Mtrix( i ) = Ansat ;
    if first_P == 0 && all( Xxyz( i,: ) ) == 1  %第一次定位成功後,紀錄當下時刻
        first_P = i ;
    end
    
end

close( bar ) ;
run_delta = run_ADR-run_trueRange;
switch obs_flag
    case { 0 }
        obs = 'Dual_EPR_LS' ;
    case { 1 }
        obs = 'Dual_PR_LS' ;
    case { 2 }
        obs = 'Dual_DSC_LS' ;
    case { 3 }
        obs = 'Dual_CSC_LS' ;
end

%---------------------資料取樣時間----------------------
sampling_time_end = runtime ;
effTimeArray = 2 : sampling_time_end ;
%---------------------顯示參數設定---------------------
[ Data_start_end_obs ] = { Data sampling_time_end obs  }

%--------------------將xyz轉成enu------------------------
rec_pos_test = mean( Xxyz(effTimeArray,:),1 ) ;      %計算STD
abs = [40 116 100];
rec_pos_act = lla2ecef( abs , 'WGS84' );
%rec_pos_act = lla2ecef([22.4259617417365 , 114.211355684335 ,53.2231195336208], 'WGS84' )   %華大測站參考位置
%rec_pos_act = [-2144855.42, 4397605.31, 4078049.85];
Xenu = zeros( runtime , 3 ) ;
NEU = xyz2llh(rec_pos_test)
for i = 1 : runtime
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    sum0 = sum0 + Xenu( i , : ).^2 ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;
std = sqrt( sum0/runtime ) ;

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
title( 'PDOP (GPS-t)' )
print( '-dpng',  'PDOP_GPS-t' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 2 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS-t)' )
print( '-dpng',  'Number of visible satellites_GPS-t' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 3 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
%axis( [-3.5 , 1.5 ,-1.5 , 2.5 ] ) ;
title( 'Scatter Plot (GPS-t)' ) ;
print( '-dpng',  'Scatter plot_GPS-t' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS-t)');
%axis( [xlim,-4,2] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-2,4] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-5,5] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS-t' , '-r600' ) ;       %Change "-r600" to the required DPI

%{
figure( 7 ) ;
plot( effTimeArray , run_delta( 14 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP14' )
print( '-dpng',  'GP14-error' , '-r600'  ) ;
figure( 8 ) ;
plot( effTimeArray , run_delta( 10 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP10' )
print( '-dpng',  'GP10-error' , '-r600'  ) ;
figure( 9 ) ;
plot( effTimeArray , run_delta( 7 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP7' )
print( '-dpng',  'GP7-error' , '-r600'  ) ;
%}
figure( 10 ) ;
plot( effTimeArray , run_recbias( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'receiver bias' )
print( '-dpng',  'receiver bias' , '-r600'  ) ;
toc