tic
format long g
clc ;
clear ;
close all ;

addpath Functions

c = 299792458 ;                                 % 光速(m/s)
L1_freq = 1575.42*10^6 ;                % L1 band frequency(Hz)
L1_lambda = c / L1_freq ;

Data = 11153 ;
obs_flag = 0 ;                                             % 0=epr ; 1=pr ; 2=dsc;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
G_B_condition = 1.05 ;
interpolation_order = 13 ;                      %內插階數

Data_List           %   從Data_List.m內讀取資料

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;
[ constellations_BEIDOU ] = goGNSS.initConstellation( 0 , 0 , 0 , 1 , 0 , 0 ) ;
%   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag,BDS_flag, QZS_flag, SBS_flag);

%-----------讀取接收機資訊------------------------
[ week_num , hwtime , useless , sat_sys , prn , cnr0 , TOW0 , Elapse_Epoch , Elapse_Code , pr , pr_rate , ADR , epr , sat_EL , sat_AZ , time ] = ...
    textread( obs_filename , '%f %f %c %c %f %f %f %f %f %f %f %f %f %f %f %f' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

%-----------initialize--------------------------
runtime = round( max( time ) ) - round ( min( time ) ) ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
EDOP_G=EDOP;                        NDOP_G=EDOP;     VDOP_G=EDOP;     HDOP_G=EDOP;     PDOP_G=EDOP;     GDOP_G=EDOP;
EDOP_B=EDOP;                         NDOP_B=EDOP;     VDOP_B=EDOP;     HDOP_B=EDOP;      PDOP_B=EDOP;     GDOP_B =EDOP;
nGsat_Mtrix = EDOP ;              nBsat_Mtrix = EDOP ;                              nsat_Mtrix = EDOP ;
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
wl_count = 0 ;     %計算GPS/BDS權重疊代次數
chg = 1 ;     % record the change of satellite amounts
rec_bias_satpos = 0 ;
Residue_set = zeros( runtime,1 ) ;
%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8

%--------------------------------讀取廣播星歷資料----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage
[ Eph_tatol_B , iono_para_B ] = RINEX_get_nav( nav_filename , constellations_BEIDOU ) ;   % load navigation massage

%----------------------------得到SP3檔案的資料---------------------------
[ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;

%---------------從ionex檔案中得到VTEC與lat,lon,time之資訊----------
[ VTEC_MAP , time_total , lat_total , lon_total ] = ...    % 注意VTEC的單位為0.1TECU
    GPS_ParseIONEXTEC( GIM_filename ) ;

%---------------solid earth tide error----------------------------------
[ time_solid , n_solid , e_solid , u_solid ] = textread( solid_filename ,'%f %f %f %f','delimiter',' ','headerlines',2); %headerlines:跳過前兩行

%----------------猜測初始位置--------------------------
llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
%llh_0 = [ 25.034214 121.564488 3 ] ;       %台北101
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
rec_bias_G = 0 ;
rec_bias_B = 0 ;

T = 1; % positioning interval
% State vector is as [x Vx y Vy z Vz b d].', i.e. the coordinate (x,y,z),
% the clock bias b, and their derivatives.

% Set f
f = @(X) ConstVelocity( X , T ) ;

% Kalman Parameters initial setting
sigma= 0.05 ;      Rhoerror = 49 ;
%Sb = (1.1e-19)*c^2 ;        Sd = (4.3e-20)*c^2 ;
%Sb = (4e-19)*c^2 ;        Sd = ((pi^2)*16e-20)*c^2 ;%The study of GPS Time Transfer Based on Extended Kalman Filter
Sb = (2e-19)*c^2 ;          Sd = ((pi^2)*56e-20)*c^2 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)
%Sb = (4e-19)*c^2 ;        Sd = (1.58e-18)*c^2 ;
X = zeros( 10 , 1 ) ;
P = eye( 10 )*5 ;

Qc = [Sb*T+Sd*T*T*T/3     Sd*T*T/2 ;
    Sd*T*T/2                Sd*T ] ;
Qxyz = sigma^2 * [T^3/3      T^2/2 ;
    T^2/2        T      ] ;
Q = blkdiag( Qxyz , Qxyz , Qxyz , Qc , Qc ) ;

first_P = 0 ;
KF = 0 ;
bar = waitbar(0,'Please wait...');
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
    
    GP_id = [] ;                        BD_id = [] ;
    Rho0_G = [] ;                    Rho0_B = [] ;
    Rho0_Gpr = [] ;                Rho0_Bpr = [] ;
    CNR_G = [] ;                     CNR_B = [] ;
    EL_G = [] ;                         EL_B = [] ;
    AZ_G = [] ;                        AZ_B = [] ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    k = 1 ;                           %每一秒都把各顆BDS的資料依序抓進來
    
    while judge == 1
        
        while sat_EL( count ) < smaller_than_elevation
            count = count + 1 ;
        end
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
            if char( sat_sys( count ) ) == 80 &&  prn( count ) ~= 193 && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999
                GP_id( j , 1) = prn( count ) ;
                CNR_G( j , 1) = cnr0( count ) ;
                Rho0_G( j , 1 ) = epr( count ) ;
                Rho0_Gpr( j , 1 ) = pr( count ) ;
                EL_G( j , 1 ) = sat_EL( count ) ;
                AZ_G( j , 1 ) = sat_AZ( count ) ;
                
                j = j+1 ;
                
            elseif char( sat_sys( count ) ) == 66 &&  prn( count ) ~= 5
                BD_id( k , 1 ) = prn( count ) ;
                CNR_B( k , 1 ) = cnr0( count ) ;
                Rho0_B( k , 1 ) = epr( count ) ;
                Rho0_Bpr( k , 1 ) = pr( count ) ;
                EL_B( k , 1 ) = sat_EL( count ) ;
                AZ_B( k , 1 ) = sat_AZ( count ) ;
                
                k = k + 1 ;
                
            else
                
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
    
    %----------------是否使用original pr----------
    if obs_flag == 1
        Rho0_G = Rho0_Gpr ;
        Rho0_B = Rho0_Bpr ;
    end
    
    %----------------每次(每一秒)定位衛星的顆數---------------
    nGsat = length( GP_id ) ;
    nBsat = length( BD_id ) ;
    nsat = nGsat + nBsat ;
    nGsat_Mtrix( i ) = nGsat ;
    nBsat_Mtrix( i ) = nBsat ;
    nsat_Mtrix( i ) = nsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nGsat ~= 0 && nBsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        %-------------------計算covariance----------------------------
        EL = [ EL_G ; EL_B ] ;
        CNR = [ CNR_G ; CNR_B ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        %-------------------選擇權重矩陣------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        
        Rr_G = Rho0_G ;
        Rr_B = Rho0_B ;
        Delta_Mtrix( 1:5 ) = 1 ;
        Residue( [1,3,5] ) = 1 ;
        Helmert = 1 ;
        
        while norm( Delta_Mtrix(1:3) ) > 1e-4 || norm( Residue([1,3,5]) ) > 0.9
            
            llh_0 = xyz2llh( rec_pos_0 ) ;
            lat_rec = llh_0( 1 ) ;          % unit : degree
            lon_rec = llh_0( 2 ) ;         % unit : degree
            height_rec = llh_0( 3 ) ;   % unit : meter
            
            Tg = zeros( nGsat , 1 ) ;              ig = zeros( nGsat , 1 ) ;
            Tb = zeros( nBsat , 1 ) ;              ib =  zeros( nBsat , 1 ) ;
            %--------------------得到transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
            [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
                satellite_positions( r_gpst , Rr_B , BD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
            
            %------------利用 UNB3 模型計算對流層誤差-------------------------
            for n = 1 : nGsat
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL_G(n) * pi/180 ) ;
                Tg( n , 1 ) = RTROP ;                   % lat_rec : radians
            end
            for n = 1 : nBsat
                [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL_B(n) * pi/180 ) ;
                Tb( n , 1 ) = RTROP_B ;
            end
            
            %------------使用GIM計算電離層誤差--------------------------
            for n = 1 : nGsat
                [ iono , VTEC ]=...
                    GIM_corr( AZ_G(n) , EL_G(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                ig( n , 1 ) = iono ;
            end
            for n = 1 : nBsat
                [ iono_B , VTEC_B ]=...
                    GIM_corr( AZ_B(n) , EL_B(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                ib( n , 1 ) = iono_B ;
            end
            
            %------------使用 time_tx 計算精密星歷位置-----------------
            %--------------initialize----------------
            sat_G = zeros( nGsat , 3 ) ;
            SCBg = zeros( nGsat , 1 ) ;
            
            %--------------------計算精密星歷衛星位置及鐘差---------------------------
            for k = 1 : nGsat
                sat_x_temp = LagrangeInter( Data_time , AX( GP_id(k),: ) , time_tx( k ) ) ;
                sat_y_temp = LagrangeInter( Data_time , AY( GP_id(k),: ) , time_tx( k ) ) ;
                sat_z_temp = LagrangeInter( Data_time , AZ( GP_id(k),: ) , time_tx( k ) ) ;
                sat_G( k , 1 ) = sat_x_temp ;
                sat_G( k , 2 ) = sat_y_temp ;
                sat_G( k , 3 ) = sat_z_temp ;
                for t = 1 : length( AX( GP_id(k),: ) )
                    if Data_time( t ) - time_tx( k ) > 0
                        Interval = t-1 : t ;
                        break
                    end
                end
                ephemeris_SCB = AS( GP_id(k),: ) ; %SCB=satellite clock bias
                SCB_temp = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
                SCBg( k , 1 ) = SCB_temp ;
            end
            
            %-------------------使用stallie_position------------------------
            sat_B = XS_tx_B ;
            SCBb = c .* dtS_B ;
            
            %------------------計算GPS相對論誤差-------------------------------
            relative_IGS = zeros( nGsat , 1 ) ;
            tgd = 0 ;
            
            for n = 1 : nGsat
                icol = find_eph( Eph_tatol , GP_id( n ) , r_gpst ) ;
                Eph = Eph_tatol( : , icol ) ;
                
                [ satp , satv ] = satellite_orbits( time_tx( n ) , Eph , icol , [] ) ;
                
                Sp3_corr = -2*dot( satp , satv )/c ;
                relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
            end
            
            %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------
            % 北斗的SCBb包含sat_bias,relative與group_delay
            Rhoc_G = Rho0_G - ig - Tg + SCBg + relative_IGS ;
            Rhoc_B = Rho0_B - ib - Tb + SCBb ;
            Rhoc = [ Rhoc_G ; Rhoc_B ] ;
            
            %---------------------------修正地球自轉--------------------------------
            [ Xs_G , Rr_G ] = Fix_Earth_Rotation( sat_G , rec_pos_0 ) ;
            [ Xs_B , Rr_B ] = Fix_Earth_Rotation( sat_B , rec_pos_0 ) ;
            Xs = [ Xs_G ; Xs_B ] ;
            
            %--------------------Extended Kalman Filter------------------------------
            if KF == 1
                g = @(X) Linearize_Measurement( X , Xs , nGsat , nBsat ) ;
                R = eye( size( Xs , 1 ) ) * Rhoerror;
                Z = Rhoc ;
                
                X0 = X ;
                [ X , P ] = Extended_Kalman_Filter( f , g , Q , R , Z , X , P ) ;
                X1 = X ;
                Residue = X1-X0 ;
                Residue_set(i) = norm( Residue([1,3,5]) ) ;
                
                rec_pos_0 = X( [1,3,5] )' ;
                Xxyz( i , : ) = rec_pos_0 ;
            end
            
            if nsat > 4 && first_P == 0
                
                [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                
                while  ( abs( V_hyber(1)/V_hyber(2) )  > G_B_condition ||  abs( V_hyber(2)/V_hyber(1) ) > G_B_condition ) && Helmert == 1
                    W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat , 1:nGsat) , ( 1/V_hyber(2) ) * W(nGsat+1:nsat , nGsat+1:nsat) ) ;
                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                    wl_count = wl_count + 1 ;
                end
                rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                rec_bias_G = rec_bias_G + Delta_Mtrix(4) ;
                rec_bias_B = rec_bias_B + Delta_Mtrix(5) ;
                rec_pos_0 = rec_pos ;
                Helmert = 0 ;
                Residue([1,3,5]) = 0 ;
            else
                Delta_Mtrix(1:3) = 0 ;
            end
            
        end
        
        if first_P == 0 && nsat < 5
            Xxyz( i , : ) = rec_pos_0 ;
        elseif first_P == 0
            Xxyz( i , : ) = rec_pos ;
            X( [1,3,5] ) = rec_pos' ;
            X(7) = rec_bias_G ;
            X(9) = rec_bias_B ;
            first_P = 1 ;
            KF = 1 ;
        end
        
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
        rec_pos_0 = Xxyz( i , : ) ;             %   將所得之位置更新為下一時刻之初始位置
    end
    %----------------時間+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    %---------------------計算固體潮汐誤差(solid tide)----------------
    solid_tide_error( i , 1 ) = interp1( time_solid , e_solid , time_day_UTC ) ;
    solid_tide_error( i , 2 ) = interp1( time_solid , n_solid , time_day_UTC ) ;
    solid_tide_error( i , 3 ) = interp1( time_solid , u_solid , time_day_UTC ) ;
    
    if i > 1 && nsat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    
    nsat_prior = nsat ;     %   紀錄上一秒衛星顆數
    
end
close( bar ) ;
switch obs_flag
    case { 0 }
        obs = 'EPR' ;
    case { 1 }
        obs = 'PR' ;
end
%---------------------資料取樣時間----------------------
effTimeArray = sampling_time_start : sampling_time_end ;

%---------------------參數設定顯示----------------------
[ Data_start_end_obs_elevation_GBcond_Q_R ] = { Data sampling_time_start sampling_time_end obs ...
    smaller_than_elevation G_B_condition Sb Sd sigma Rhoerror }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;

Xenu = zeros( runtime , 3 ) ;
for i = 1 : runtime
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( var( Xenu(effTimeArray,:))  ) ;
rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
STD = [ std_enu' ; std_en' ] ;

rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ;
RMS = [ rms_enu' ; rms_en' ] ;

error = [ std' , rms' ]
Error = [ STD , RMS ]
%{
%--------------------GPS/BDS DOP--------------------------
    DOP_mean = [ mean( HDOP ) mean( VDOP ) mean( PDOP ) ] ;
%{
    figure( 1 ) ;
    plot( effTimeArray , HDOP( effTimeArray ) , 'r-' , ...
            effTimeArray , VDOP( effTimeArray ) , 'b-' , ...
            effTimeArray , PDOP( effTimeArray ) , 'g-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;                                   ylabel( 'DOP value' ) ;
    legend( 'HDOP' , 'VDOP' , 'PDOP' ) ;    title( 'GPS/BDS DOP' )
%}
    
    figure( 2 ) ;
    plot( effTimeArray , GDOP( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'GDOP value' ) ;
    legend( 'GPS/BDS' ) ;
    title( 'GPS/BDS GDOP' )
    print( '-dpng',  'GDOP_Daul' , '-r600'  ) ;
    
    figure( 3 ) ;
    plot( effTimeArray , GDOP_G( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'GDOP value' ) ;
    legend( 'GPS' ) ;
    title( 'GPS GDOP' )
    print( '-dpng',  'GDOP_GPS' , '-r600'  ) ;
    
    figure( 4 ) ;
    plot( effTimeArray , GDOP_B( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'GDOP value' ) ;
    legend( 'GPS' ) ;
    title( 'BDS GDOP' )
    print( '-dpng',  'GDOP_BDS' , '-r600'  ) ;
        
    figure( 5 ) ;
    plot( effTimeArray , GDOP( effTimeArray ) , 'r-' ,...
        effTimeArray , GDOP_G( effTimeArray ) , 'b-' , ...
        effTimeArray , GDOP_B( effTimeArray ) , 'g-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'GDOP value' ) ;
    legend( 'GPS/BDS' , 'GPS' , 'BDS' ) ;
    title( 'GDOP' )
    print( '-dpng',  'GDOP_three systems' , '-r600'  ) ;
    
%--------------------Number of visible satellites----------------------
    figure( 6 ) ;
    plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'number' ) ;
    ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
    legend( 'Number of visible satellites' ) ;
    title( 'GPS/BDS satellite' )
    print( '-dpng',  'Number of Dual' , '-r600'  ) ;

    figure( 7 ) ;
    plot( effTimeArray , nGsat_Mtrix( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'number' ) ;
    legend( 'Number of visible satellites' ) ;
    title( 'GPS satellte' ) ;
    print( '-dpng',  'Number of GPS' , '-r600'  ) ;

    figure( 8 ) ;
    plot( effTimeArray , nBsat_Mtrix( effTimeArray ) , 'b-' ) ;
    grid on ;
    xlabel( 'time(s)' ) ;
    ylabel( 'number' ) ;
    legend( 'Number of visible satellites' ) ;
    title( 'BDS satellite' ) ;
    print( '-dpng',  'Number of BDS' , '-r600'  ) ;
%}
figure( 1 ) ;
plot( effTimeArray , Residue_set(effTimeArray) , 'b' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Residue (m)') ;
title( '量測餘量' ) ;
print( '-dpng',  'KF量測餘量' , '-r600'  ) ;       %Change "-r600" to the required DPI

%--------------------2D plot----------------------
figure( 10 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
title( 'Scatter Plot (GPS/BDS)' ) ;
%axis( [ -7,16,-10,10 ] );
print( '-dpng',  'Scatter plot_KF' , '-r600'  ) ;       %Change "-r600" to the required DPI
%{
%-------------------3D plot------------------------------
    figure( 11 );
    plot(effTimeArray , Xenu( effTimeArray , 1 ) , 'r-.' ,...
            effTimeArray , Xenu( effTimeArray , 2 ) , 'b-.' ,...
            effTimeArray , Xenu( effTimeArray , 3 ) , 'g-.' ) ;
    grid on;
    xlabel('time(s)');
    ylabel('Errors(m)');
    legend('East','North','UP');
    print( '-dpng',  'ENU error' , '-r600'  ) ;
%}
figure( 12 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS/BDS)');
%axis( [xlim,-2,4] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-10,0] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-20,20] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_KF' , '-r600' ) ;       %Change "-r600" to the required DPI
toc