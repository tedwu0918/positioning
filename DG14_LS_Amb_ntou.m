format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % 光速(m/s)
lambda = c /(1575.42*10^6) ;

Data = 0427 ;
obs_flag = 1 ;
converge_time = 3000000 ;
initial_time = converge_time + 20 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
interpolation_order = 13 ;                      %內插階數

switch Data
    case{ 0427 }
        addpath DG14_0427
        year = 2016 ;     week_num = 1894 ;         month = 4 ;      day = 27 ;      day_of_year = 118 ;
        obs_filename = 'NTOU0427test.16o' ;                                 % DG14的數據檔案
        SP3_filename = 'igu18943_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg1180.16i' ;                                   % GIM全球電離層 ionex檔名                         
        nav_filename = 'brdm1180.16p' ;                                  % 廣播星歷導航資料
end

%-----------讀取接收機資訊------------------------
[ TOW , Prn , L1 , C1 ] = readOBS( obs_filename ) ;
 
%-----------initialize--------------------------
SVMax = 32 ;
runtime = TOW(end)-TOW(1) ;
%runtime = 1000 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
wl_count = 0 ;     %計算GPS/BDS權重疊代次數
chg = 1 ; % record the change of satellite amounts
PRr = zeros( 32,2 ) ; 
Carrier = zeros( 32,2 ) ;
Smooth_Count = zeros( 32,1 ) ;
Smooth_Code = zeros( 32,1 ) ; 
nsat_Mtrix = zeros( runtime,1 ) ; 
N_set = zeros( runtime,1 ) ;

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;

%--------------------------------讀取廣播星歷資料----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage

%--------------------------------得到SP3檔案的資料----------------------------------
[ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
    SP3_read( year , month , day , SP3_filename , week_num , TOW(1) , interpolation_order ) ;

%---------------從ionex檔案中得到VTEC與lat,lon,time之資訊----------
[ VTEC_MAP , time_total , lat_total , lon_total ] = ...    % 注意VTEC的單位為0.1TECU
            GPS_ParseIONEXTEC( GIM_filename ) ;

%----------------猜測初始位置--------------------------
rec_pos_0 = [ -3042325.7824 , 4911006.2700 , 2694068.1809  ] ;
rec_bias = 0 ;
bar = waitbar(0,'Please wait...') ;
count = 1 ;

time_week = TOW(1) - 1 ;

for i = 1 : runtime
    
    time_week = time_week +1 ;
    
    str=['Positinoing ',num2str(i),'s'];    
    waitbar( i/runtime , bar , str )   % computation here %

    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
                    SP3_read( year , month , day , SP3_filename , week_num , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
    
    rec_bias_satpos = 0 ;  
    id = [] ;    
    Phi0_temp= [] ;
    Rho0_temp = [] ;  
    id_temp = [] ;
    EL = [] ;        
    AZI = [] ;           
    PR = zeros( 32,1 ) ;
    Carrier( :,1 ) = Carrier( :,2 ) ;
    Carrier( :,2 ) = 0 ;
    SVxyz = zeros( 1,3 ) ;
    %----------------------------------------------- 
    j = 1 ;
    while TOW( count ) == time_week ;        
        if all( AX( Prn(count),: ) ) == 1 && AS( Prn(count),1 ) ~= 999999.999999
            id_temp(j,1) = Prn( count ) ;
            Phi0_temp(j,1) = L1( count ) * lambda ;
            Rho0_temp(j,1) = C1( count ) ;
            j = j + 1 ;
        end
        count = count+1 ;
    end
       
    if size(id_temp,1) == 0
        Xxyz( i,: ) = Xxyz( i-1,: ) ;
        continue
    end
    
    %-----------得到接收時間---------------------
    r_gpst = weektow2time ( week_num , time_week , 'G' ) ;  %GPST
    
    llh_0 = xyz2llh( rec_pos_0 ) ;
    lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
    lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
    height_rec = llh_0( 3 ) ;                  % unit : meter
       
    nsat = size(id_temp,1) ;
    Tropo = zeros( nsat , 1 ) ;
    Iono = zeros( nsat , 1 ) ;
    %--------------------得到transmission time ( time_tx )----------------------
    [ XS , dtS , XS_tx , VS_tx , time_tx_temp , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
        satellite_positions( r_gpst , Rho0_temp , id_temp , Eph_tatol , [] , [] , Tropo , Iono , rec_bias_satpos ) ;    %注意,迭代時要帶入地球修正後的虛擬距離(Rr)
    
    %--------------------計算精密星歷衛星位置及鐘差---------------------------
    sat_temp = zeros( nsat,3 ) ;
    SCBg_temp = zeros( nsat,1 ) ;
    for k = 1 : nsat
        sat_temp( k , 1 ) = LagrangeInter( Data_time , AX( id_temp(k),: ) , time_tx_temp( k ) ) ;
        sat_temp( k , 2 ) = LagrangeInter( Data_time , AY( id_temp(k),: ) , time_tx_temp( k ) ) ;
        sat_temp( k , 3 ) = LagrangeInter( Data_time , AZ( id_temp(k),: ) , time_tx_temp( k ) ) ;
        for t = 1 : length( AX( id_temp(k),: ) )
            if Data_time( t ) - time_tx_temp( k ) > 0
                Interval = t-1 : t ;
                break
            end
        end
        ephemeris_SCB = AS( id_temp(k),: ) ; %SCB=satellite clock bias
        SCBg_temp( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx_temp( k ) ) * c * 1e-6 ;
    end
    
    EL_temp = zeros( nsat,1) ;
    AZI_temp = zeros( nsat,1 ) ;
    for k = 1 : nsat
        SVxyz = sat_temp(k,:) ;
        %----------------計算仰角&方位角,並轉成角度單位---------------
        [ EL_temp(k,1) , AZI_temp(k,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , SVxyz ) ;
    end
    EL_temp = EL_temp * ( 180/pi ) ;
    AZI_temp = AZI_temp * ( 180/pi ) ;
    
    
    clear Phi0 Rho0 id sat time_tx EL AZI
    j = 1 ;
    for k = 1 :nsat
        if EL_temp(k) >= smaller_than_elevation
            Phi0(j,1) = Phi0_temp(k) ;
            Rho0(j,1) = Rho0_temp(k) ;
            id(j,1) = id_temp(k) ;
            sat( j,: ) = sat_temp( k,: ) ;
            time_tx(j,1) = time_tx_temp(k) ;
            EL(j,1) = EL_temp(k) ;
            AZI(j,1) = AZI_temp(k) ;
            j = j + 1 ;
        end
    end
    
    for k = 1 : size(id,1)
        PR( id(k),1 ) = Rho0(k) ;
        Carrier( id(k),1 ) =Phi0(k) ;    
    end
    
    %----------------每次(每一秒)定位衛星的顆數---------------
    nsat = length( id ) ;
    nsat_Mtrix( i ) = nsat ;
    
    %-----------------------------------------------------------------------------------------------    
    if nsat ~= 0                       
        %----------------選擇觀測量------------------
        if obs_flag == 3
            [ Rho0 , Smooth_Code , Smooth_Count ] = Carrier_Smoothed_Code( Carrier , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        end
        
        Rr = Rho0 ;
        Delta_Mtrix(1:3) = 1 ;
        W = eye( nsat ) ;
        
        while norm( Delta_Mtrix(1:3) ) > 0.01
            
            EL = zeros( nsat,1) ;
            AZI = zeros( nsat,1 ) ;
            for k = 1 : nsat
                SVxyz = sat(k,:) ;
                %----------------計算仰角&方位角,並轉成角度單位---------------
                [ EL(k,1) , AZI(k,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , SVxyz ) ;
            end
            EL = EL * ( 180/pi ) ;
            AZI = AZI * ( 180/pi ) ;            
            
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
            Phic = Phi0 + Iono - Tropo +SCBg + relative_IGS ;
            
            %---------------------------修正地球自轉--------------------------------
            [ Xs , Rr ] = Fix_Earth_Rotation( sat , rec_pos_0 ) ;
            
            %-----------------------最小平方法定位-------------------------   
            if nsat > 3     %   單系統至少需要4顆衛星來定位
                [ Delta_Mtrix , Delta_Rho , Delta_Phi , Nf , G ] = LS_Amb( Xs , rec_pos_0 , Rhoc , Phic , nsat , rec_bias ) ;
                
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
        end
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
        rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置        
        
        for k = 1 : size(id,1)
            if id(k) == 20
                N_set( i ) = Nf(k) ;
            end
        end
        
    end
    
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
effTimeArray = 10 : sampling_time_end ;

%---------------------顯示參數設定---------------------
%[ Data_start_end_obs_elevation_GBcond ] = { Data sampling_time_start sampling_time_end obs  }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;      %計算STD

Xenu = zeros( runtime , 3 ) ;
for i = 1 : runtime 
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( var( Xenu(effTimeArray,:))  ) ;
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

figure( 3 ) ;
plot( effTimeArray , PDOP( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS)' )
print( '-dpng',  'PDOP_GPS' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 6 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS)' )
print( '-dpng',  'Number of visible satellites_GPS' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 10 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 12 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI
