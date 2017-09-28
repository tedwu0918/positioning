tic
format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % 光速(m/s)
lambda = c /(1575.42*10^6) ;

Data = 0118 ;
obs_flag = 3 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 30 ;
chenckcount = 30 ;
initial_time = converge_time + 10 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
interpolation_order = 13 ;                      %內插階數

switch Data
    case{ 0109 }
        addpath DG14_0109
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1931 ;         month = 1 ;      day = 9 ;      day_of_year = 9 ;
        obs_filename = '_msvgpsuv_com4_0109.csv' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19311_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0090.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0090.17p' ;                                  % 廣播星歷導航資料
        TDD_filename = 'bout.txt'                                           %tripple difference data
    case{ 0116 }
        addpath DG14_0116
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1932 ;         month = 1 ;      day = 16 ;      day_of_year = 16 ;
        obs_filename = '1-16-s.txt' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19321_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0160.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0160.17p' ;                                  % 廣播星歷導航資料
        TDD_filename = 'bout.txt'     ;                                      %tripple difference data
    case{ 0118 }
        addpath DG14_0118
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1932 ;         month = 1 ;      day = 18 ;      day_of_year = 18 ;
        obs_filename = '1-18-s2.txt' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19323_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0180.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0180.17p' ;                                  % 廣播星歷導航資料
        TDD_filename = 'bout182.txt'     ;                                      %tripple difference data
    case{ 0208 }
        addpath DG14_0208
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1935 ;         month = 2 ;      day = 8 ;      day_of_year = 39 ;
        obs_filename = '2-8-s.txt' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19353_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0390.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0390.17p' ;                                  % 廣播星歷導航資料
        TDD_filename = 'bout.txt'     ;                                      %tripple difference data
    case{ 0209 }
        addpath DG14_0208
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1935 ;         month = 2 ;      day = 9 ;      day_of_year = 40 ;
        obs_filename = '2-9-s.txt' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19353_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0390.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0390.17p' ;                                  % 廣播星歷導航資料
        TDD_filename = 'bout.txt'     ;                                      %tripple difference data
        
end

%-----------讀取接收機資訊------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

%-----------initialize--------------------------
SVMax = 32 ;
runtime = round( max( time ) ) - round ( min( time ) ) +1;
runtime = 180 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
nsat_Mtrix = EDOP ;
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
Residue_set = zeros( runtime,1 ) ;
wl_count = 0 ;     %計算GPS/BDS權重疊代次數
chg = 1 ; % record the change of satellite amounts
PRr = zeros( 32,2 ) ;
Carrier = zeros( 32,2 ) ;
Smooth_Count = zeros( 32,1 ) ;
Smooth_Code = zeros( 32,1 ) ;
first_postime = 0 ;
checkmatrix = zeros( 32,3 ) ;

%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;

%--------------------------------讀取廣播星歷資料----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage

%----------------------------得到SP3檔案的資料---------------------------
[ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;

%---------------從ionex檔案中得到VTEC與lat,lon,time之資訊----------
[ VTEC_MAP , time_total , lat_total , lon_total ] = ...    % 注意VTEC的單位為0.1TECU
    GPS_ParseIONEXTEC( GIM_filename ) ;

%----------------猜測初始位置--------------------------
llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
rec_bias = 0 ;

%----------------KF setting--------------------------
T = 1; % positioning interval
% State vector is as [x Vx y Vy z Vz b d].', i.e. the coordinate (x,y,z),
% the clock bias b, and their derivatives.


% Set f
f = @(X) ConstVelocity( X , T ) ;

% Kalman Parameters initial setting
sigma= 1 ;      Rhoerror = 5 ;       Pvalue = 10 ;
%Sb = (1.1e-19)*c^2 ;        Sd = (4.3e-20)*c^2 ;
%Sb = (4e-19)*c^2 ;        Sd = ((pi^2)*16e-20)*c^2 ;%The study of GPS Time Transfer Based on Extended Kalman Filter
Sb = (1.1e-19)*c^2 ;          Sd = ((pi^2)*56e-21)*c^2 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)
%Sb = (4e-19)*c^2 ;        Sd = (1.58e-18)*c^2 ;

X = zeros( 8 , 1 ) ;
P = eye( 8 )*Pvalue ;
Xu = zeros( 8 , 1 ) ;

Qc = [Sb*T+Sd*T*T*T/3     Sd*T*T/2 ;
    Sd*T*T/2                Sd*T ] ;
Qxyz = sigma * [T^3/3      T^2/2 ;
    T^2/2        T      ] ;
Q = blkdiag( Qxyz , Qxyz , Qxyz , Qc ) ;

first_P = 1 ;
KF = 0 ;
%-------------------------------------------------------------------------------------
bar = waitbar(0,'Please wait...') ;

for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
    
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
            if prn( count ) ~= 193 && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999 %&&prn( count ) ~= 6 &&prn( count ) ~= 24 
                
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
                    checkmatrix( id( j , 1),1 ) = 1 ;
                    checkmatrix( id( j , 1),3 ) = checkmatrix( id( j , 1),3 ) + 1 ;
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
        if obs_flag == 3 && chenckcount ~=0
            [W , checkmatrix ] = adjweightingfunc( W , checkmatrix , id , chenckcount ) ;
        end
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
            %--------------------得到transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                satellite_positions( r_gpst , Rr , id , Eph_tatol , [] , [] , Tropo , Iono , rec_bias_satpos ) ;    %注意,迭代時要帶入地球修正後的虛擬距離(Rr)
            
            %------------利用 UNB3 模型計算對流層誤差-------------------------
            for n = 1 : nsat
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL(n)*pi/180 ) ;
                Tropo( n , 1 ) = RTROP ;
            end
            
            %------------使用GIM計算電離層誤差--------------------------
            for n = 1 : nsat
                [ iono , VTEC ]=...
                    GIM_corr( AZI(n) , EL(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                Iono( n , 1 ) = iono ;
            end
            
            %------------使用 time_tx 計算精密星歷位置-----------------
            %--------------initialize----------------
            sat = zeros( nsat , 3 ) ;
            SCBg = zeros( nsat , 1 ) ;
            %--------------------計算精密星歷衛星位置及鐘差---------------------------
            for k = 1 : nsat
                sat( k , 1 ) = LagrangeInter( Data_time , AX( id(k),: ) , time_tx( k ) ) ;
                sat( k , 2 ) = LagrangeInter( Data_time , AY( id(k),: ) , time_tx( k ) ) ;
                sat( k , 3 ) = LagrangeInter( Data_time , AZ( id(k),: ) , time_tx( k ) ) ;
                for t = 1 : length( AX( id(k),: ) )
                    if Data_time( t ) - time_tx( k ) > 0
                        Interval = t-1 : t ;
                        break
                    end
                end
                ephemeris_SCB = AS( id(k),: ) ; %SCB=satellite clock bias
                SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
            end
            
            %------------------計算GPS相對論誤差-------------------------------
            relative_IGS = zeros( nsat , 1 ) ;
            tgd = 0 ;
            
            for n = 1 : nsat
                icol = find_eph( Eph_tatol , id( n ) , r_gpst ) ;
                Eph = Eph_tatol( : , icol ) ;
                
                [ satp , satv ] = satellite_orbits( time_tx( n ) , Eph , icol , [] ) ;
                
                Sp3_corr = -2*dot( satp , satv )/c ;
                relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
            end
            
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
        
        if first_P == 1 && nsat < 4
            Xxyz( i , : ) = rec_pos_0 ;
        elseif first_P == 1
            Xxyz( i , : ) = rec_pos ;
            X( [1,3,5] ) = rec_pos' ;
            X(7) = Delta_Mtrix(4) ;
            KF = 1 ;
            first_P = 0 ;
            if first_postime ==0
                first_postime = i
            end
        end
        
        %--------------------Kalman Filter------------------------------
        
        if KF == 1 %&& i > 1 && any( Xxyz(i-1,:) ) == 1
            %for tt = 1 : 10
                g = @(X) Linearize_Measurement_single( X , Xs , nsat) ;
                R = eye( size( Xs , 1 ) ) * Rhoerror;
                Z = Rhoc ;
                
                X0 = X ;
                [ X , P ] = Extended_Kalman_Filter( f , g , Q , R , Z , X , P ) ;
                X1 = X ;
                Residue = X1-X0 ;
                Residue_set(i) = norm( Residue([1,3,5]) ) ;
                
                rec_pos_0 = X( [1,3,5] )' ;
                Xxyz( i , : ) = rec_pos_0 ;
            %end
            Xxyz( i , : ) = rec_pos ;
            
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

%----------------------處理三次差分-----------------------
[  TDx , TDy , TDz , s0 , s1 , s2 ] = ...
    textread( TDD_filename , '%f %f %f %f %f %f  '  ) ;

temp_TDxyz = [TDx  TDy  TDz] ;

%---------------------資料取樣時間----------------------
if obs_flag == 3  %g使用CSC時門檻時間前不列入計算
    temp_TDxyz(1:(initial_time - 2) , :) = [] ;
    Xxyz(1:initial_time , :) = [] ;
    nsat_Mtrix(1:initial_time ) = [] ;
    sampling_time_end = runtime - initial_time ;
    TimeArray = first_postime : sampling_time_end ;
    effTimeArray_error = TimeArray ;
else
    Xxyz(1 : 2 , :) = [] ;
    nsat_Mtrix(1 : 2) = [] ;
    sampling_time_end = runtime - 2 ;
    TimeArray = first_postime : sampling_time_end ;
    effTimeArray_error = round(length(TimeArray)/100) -1 + first_postime : sampling_time_end ;
end

Time = length(TimeArray) ;

%---------------------顯示參數設定---------------------
[ Data_start_end_obs ] = { Data sampling_time_end obs  }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray_error,:),1 ) ;      %計算STD
%rec_pos_act = lla2ecef( [25.149442 , 121.777873 , 73.6] , 'WGS84' );
Xenu = zeros( Time , 3 ) ;
MSE_enu = zeros(1,3);
trajectory = zeros( Time , 3 ) ;
for i = first_postime : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    
end
trajectory = bsxfun(@minus , Xenu , Xenu(sampling_time_end , :) ) ;
TDxyz = bsxfun(@minus , temp_TDxyz , temp_TDxyz( sampling_time_end , : ) ) ;
for i = first_postime : sampling_time_end
    error( i , : ) = trajectory( i  , : ) - TDxyz( i , : ) ;
    MSE_enu = MSE_enu + error( i , : ).^2 ;
    
end

%--------------------計算STD & RMS------------------------
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;
std = sqrt(Time \ MSE_enu) ;
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
plot( TimeArray , PDOP( TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS)' )
print( '-dpng',  'PDOP_GPS' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 2 ) ;
plot( TimeArray , nsat_Mtrix( TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(TimeArray) )-0.5  max( nsat_Mtrix(TimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS)' )
print( '-dpng',  'Number of visible satellites_GPS' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 3 ) ;
plot( Xenu( TimeArray , 1 ) , Xenu( TimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
%axis( [-4 , 4 ,-4 , 4 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( TimeArray , error( TimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(TimeArray, error( TimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(TimeArray, error( TimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-20,20] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI

figure( 5 ) ;
plot( trajectory( TimeArray , 1 ) , trajectory( TimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East (m)' ) ;
ylabel( 'North (m)') ;
axis equal ;
%axis( [-4 , 4 ,-4 , 4 ] ) ;
title( 'trajectory ' ) ;

hold on
plot( TDxyz( TimeArray , 1 ) , TDxyz( TimeArray , 2 ) , 'r.' ) ;
print( '-dpng',  'trajectory_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI


toc