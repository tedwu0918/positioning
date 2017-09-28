tic
format long g
clc ;
clear ;
close all ;

addpath Functions

c = 299792458 ;                                 % 光速(m/s)
L1_freq = 1575.42*10^6 ;                % L1 band frequency(Hz)
L1_lambda = c / L1_freq ;

Data = 2151 ;
no_solid_tide = 2151 ;
Pr_time = 600 ;
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
Xxyz = zeros( runtime , 6 ) ;
Bias_set = zeros( runtime , 3 ) ;
wl_count = 0 ;     %計算GPS/BDS權重疊代次數
chg = 1 ; % record the change of satellite amounts
rec_bias = zeros( 2,2 ) ; 
xyz = zeros( 2,3 ) ;
first = 0 ;         %   用來記錄第一次可定位的時間
SV_Center = zeros( runtime,3 ) ;

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
if Data ~= no_solid_tide
    [ time_solid , n_solid , e_solid , u_solid ] = textread( solid_filename ,'%f %f %f %f','delimiter',' ','headerlines',2); %headerlines:跳過前兩行
end
%----------------猜測初始位置--------------------------
llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
xyz( 1,: ) = rec_pos_0 ;
xyz( 2,: ) = rec_pos_0 ;
SV_ref = 1 ;

bar = waitbar(0,'Please wait...') ;
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];    
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
                    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
         
    rec_bias_satpos = 0 ;
    
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
       
    %----------------每次(每一秒)定位衛星的顆數---------------
    nGsat = length( GP_id ) ;    
    nBsat = length( BD_id ) ;
    nsat = nGsat + nBsat ;
    nGsat_Mtrix( i ) = nGsat ;
    nBsat_Mtrix( i ) = nBsat ;    
    nsat_Mtrix( i ) = nsat ; 
    
    %----------------先使用EPR定位,後使用PR定位---------------
    for y = 1:2     % 1:EPR 2:PR
       
        if y == 1            
            rec_pos_0 = xyz( y,: ) ;
        elseif y == 2 && i <= Pr_time           
            rec_pos_0 = xyz( y,: ) ;
            Rho0_G = Rho0_Gpr ;
            Rho0_B = Rho0_Bpr ;
        elseif y == 2 && i > Pr_time
            break
        end
        
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
            %W = eye( nsat ) ;
            
            Rr_G = Rho0_G ;
            Rr_B = Rho0_B ;
            Delta_Mtrix( 1:5 ) = 1 ;
            Helmert = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 1e-2
                
                llh_0 = xyz2llh( rec_pos_0 ) ;
                lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
                lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
                height_rec = llh_0( 3 ) ;                  % unit : meter
                
                Tg = zeros( nGsat , 1 ) ;
                ig = zeros( nGsat , 1 ) ;
                Tb = zeros( nBsat , 1 ) ;
                ib =  zeros( nBsat , 1 ) ;
                
                %--------------------得到transmission time ( time_tx )----------------------
                [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                    satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
                [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
                    satellite_positions( r_gpst , Rr_B , BD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
                
                %------------利用 UNB3 模型計算對流層誤差-------------------------
                for n = 1 : nGsat
                    [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL_G(n)*pi/180 ) ;
                    Tg( n , 1 ) = RTROP ;
                end
                for n = 1 : nBsat
                    [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL_B(n)*pi/180 ) ;
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
                    sat_G( k , 1 ) = LagrangeInter( Data_time , AX( GP_id(k),: ) , time_tx( k ) ) ;
                    sat_G( k , 2 ) = LagrangeInter( Data_time , AY( GP_id(k),: ) , time_tx( k ) ) ;
                    sat_G( k , 3 ) = LagrangeInter( Data_time , AZ( GP_id(k),: ) , time_tx( k ) ) ;
                    for t = 1 : length( AX( GP_id(k),: ) )
                        if Data_time( t ) - time_tx( k ) > 0
                            Interval = t-1 : t ;
                            break
                        end
                    end
                    ephemeris_SCB = AS( GP_id(k),: ) ; %SCB=satellite clock bias
                    SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
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
                
                %-----------------------最小平方法定位-------------------------
                if nsat > 4
                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias( y,1 ) , rec_bias( y,2 ) ) ;
                    
                    while  ( abs( V_hyber(1)/V_hyber(2) )  > G_B_condition ||  abs( V_hyber(2)/V_hyber(1) ) > G_B_condition ) && Helmert == 1
                        W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat , 1:nGsat) , ( 1/V_hyber(2) ) * W(nGsat+1:nsat , nGsat+1:nsat) ) ;
                        
                        [ Delta_Mtrix , Delta_Rho , V_hyber ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias( y,1 ) , rec_bias( y,2 ) ) ;
                        wl_count = wl_count + 1 ;
                    end
                    rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias( y,1 ) = rec_bias( y,1 ) + Delta_Mtrix(4) ;
                    rec_bias( y,2 ) = rec_bias( y,2 ) + Delta_Mtrix(5) ;
                    rec_pos_0 = rec_pos ;
                    Helmert = 0 ;
                else
                    break
                end
                
            end
            
            if nsat < 5
                Xxyz( i , 3*y-2:3*y ) = rec_pos_0 ;
            else
                Xxyz( i , 3*y-2:3*y ) = rec_pos ;
            end
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , 1:3 ) ) ;
            rec_pos_0 = Xxyz( i , 3*y-2:3*y ) ;            %   將所得之位置更新為下一時刻之初始位置
            xyz( y,: ) = rec_pos_0 ;
            %-----------------------計算幾何中心----------------------------
            SV_Center(i,:) = mean( Xs ) ;
            
            if first == 0   
                ft = i ;    %   紀錄第一次可定位的時間( first positioning )
                first = 1 ;
            end
        
        end        
    end
    
    if i == Pr_time         
        Pr_center = mean( Xxyz( ft:Pr_time, 4:6 ) , 1 ) ;           %   求PR中心點
        Epr_center = mean( Xxyz( ft:Pr_time, 1:3 ) , 1 ) ;         %   求EPR中心點
        
        Bias = Epr_center - Pr_center ;     %   利用EPR中心與PR中心計算得到Bias
        Xxyz( ft:Pr_time, 1:3 ) = bsxfun( @minus , Xxyz( ft:Pr_time, 1:3 ) , Bias ) ;   %   利用Bias修正EPR位置
        Bias_set( ft:Pr_time, : ) = bsxfun( @minus , Bias , zeros(  Pr_time-ft+1, 3  ) ) ;
        
    elseif i > Pr_time && nsat ~= nsat_prior
        Bias = Xxyz( i , 1:3 ) - Pr_center ;
        Xxyz( i , 1:3 ) = Xxyz( i , 1:3 ) - Bias ;
        Bias_set( i,: ) = Bias ;
    elseif  i > Pr_time && nsat == nsat_prior
        Xxyz( i , 1:3 ) = Xxyz( i , 1:3 ) - Bias ;
        Bias_set( i,: ) = Bias ;
    end
    
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
        if i > 1 && nsat ~= nsat_prior
            nsat_change(chg) = i ;
            chg = chg + 1 ;
        end
        
        if i > 1 && nsat ~= nsat_prior && SV_ref == 1
            SV_Center_Ref = SV_Center( i-1,: ) ;
            SV_ref = 0 ;
        end
            
        nsat_prior = nsat ; %   紀錄上一秒衛星顆數
        
end
close( bar ) ;
%---------------------資料取樣時間----------------------
sampling_time_end = runtime ;
effTimeArray = sampling_time_start : sampling_time_end ;

%---------------------顯示參數設定----------------------
[ Data_start_end_obs_elevation_GBcond ] = { Data sampling_time_start sampling_time_end smaller_than_elevation G_B_condition }

%--------------------將xyz轉成enu------------------------
if Data==122
    rec_pos_act = mean( Xxyz(effTimeArray,1:3),1 ) ;
end

Xenu = zeros( runtime , 6 ) ;
for y = 1:2    
    
    for i = 1 : runtime
        Xenu( i , 3*y-2:3*y ) = xyz2enu( Xxyz( i , 3*y-2:3*y  ) , rec_pos_act )' ;
    end
    Xenu( : , 3*y-2:3*y  ) = Xenu( : , 3*y-2:3*y  ) - solid_tide_error ;
    
    %--------------------計算STD & RMS------------------------
    std = sqrt( var( Xenu( effTimeArray,3*y-2:3*y ))  ) ;
    rms = bsxfun( @power , mean( Xenu( effTimeArray , 3*y-2:3*y ).^2,1 ) , 0.5 ) ;
    
    %-----------------四捨五入(round)取到小數點第X位---------------------
    std_enu = round( std .* 1e4 ) ./ 1e4 ;
    std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
    STD = [ std_enu' ; std_en' ] ;
    
    rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
    rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ;
    RMS = [ rms_enu' ; rms_en' ] ;
    
    error = [ std' , rms' ]
    Error = [ STD , RMS ]
end
%{
%--------------------Geometric center of the satellite constellation---------------------
figure( 1 ) ;
plot( effTimeArray , SV_Center( effTimeArray , 1 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'X (m)' ) ;
title( 'Geometric Center - X' )
print( '-dpng',  '幾何中心X' , '-r600' ) ;       %Change "-r600" to the required DPI

figure( 2 ) ;
plot( effTimeArray , SV_Center( effTimeArray , 2 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Geometric Center - Y' )
print( '-dpng',  '幾何中心Y' , '-r600' ) ;       %Change "-r600" to the required DPI

figure( 3 ) ;
plot( effTimeArray , SV_Center( effTimeArray , 3 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Z (m)' ) ;
title( 'Geometric Center - Z' )
print( '-dpng',  '幾何中心Z' , '-r600' ) ;       %Change "-r600" to the required DPI

figure( 4 ) ;
plot( effTimeArray , Bias_set( effTimeArray,1 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Bias (m)' ) ;
title( 'Bias of X' ) ;
print( '-dpng',  'Bias of X' , '-r600'  ) ;

figure( 5 ) ;
plot( effTimeArray , Bias_set( effTimeArray,2 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Bias (m)' ) ;
title( 'Bias of Y' ) ;
print( '-dpng',  'Bias of Y' , '-r600'  ) ;

figure( 6 ) ;
plot( effTimeArray , Bias_set( effTimeArray,3 ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Bias (m)' ) ;
title( 'Bias of Z' ) ;
print( '-dpng',  'Bias of Z' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 7 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS/BDS)' )
print( '-dpng',  'Number of Dual' , '-r600'  ) ;
%}
%{
figure( 7 ) ;
plot( effTimeArray , nGsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
title( 'Number of Visible Satellites (GPS)' ) ;
print( '-dpng',  'Number of GPS' , '-r600'  ) ;

figure( 8 ) ;
plot( effTimeArray , nBsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'number' ) ;
title( 'Number of Visible Satellites (BDS)' ) ;
print( '-dpng',  'Number of BDS' , '-r600'  ) ;

%}


figure( 9 ) ;
plot( effTimeArray , PDOP( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS/BDS)' )
print( '-dpng',  'PDOP_Dual' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 10 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East (m)' ) ;
ylabel( 'North (m)') ;
axis equal ;
title( 'Scatter Plot (GPS/BDS)' ) ;
%axis( [ -7,16,-10,10 ] );
print( '-dpng',  'Scatter plot_Dual_LS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 11 ) ;
plot( Xenu( sampling_time_start:Pr_time , 4 ) , Xenu( sampling_time_start:Pr_time , 5 ) , 'r.' ,...
    Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.') ;
grid on ;
xlabel( 'East (m)' ) ;
ylabel( 'North (m)') ;
axis equal ;
title( 'Scatter Plot (GPS/BDS)' ) ;
legend( '600-second PR','EPR' ) ;
%axis( [ -7,16,-10,10 ] );
print( '-dpng',  'Scatter plot_AEPR' , '-r600'  ) ;       %Change "-r600" to the required DPI
%{
%-------------------3D plot------------------------------
figure( 11 );
plot(effTimeArray , Xenu( effTimeArray , 1 ) , 'r-.' ,...
    effTimeArray , Xenu( effTimeArray , 2 ) , 'b-.' ,...
    effTimeArray , Xenu( effTimeArray , 3 ) , 'g-.' ) ;
grid on;
xlabel('Time (s)');
ylabel('Error (m)');
legend('East','North','UP');
print( '-dpng',  'ENU error_Dual_LS' , '-r600'  ) ;
%}
figure( 12 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS/BDS)');
%axis( [xlim,-2,4] );
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-10,0] );
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
grid on;                        hold on;
%axis( [xlim,-20,20] );
print( '-dpng',  'ENU_Dual_LS' , '-r600' ) ;       %Change "-r600" to the required DPI

toc