%-----------initialize--------------------------
GMax = 32 ;
BMax = 15 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
EDOP_G=EDOP;                        NDOP_G=EDOP;     VDOP_G=EDOP;     HDOP_G=EDOP;     PDOP_G=EDOP;     GDOP_G=EDOP;
EDOP_B=EDOP;                         NDOP_B=EDOP;     VDOP_B=EDOP;     HDOP_B=EDOP;      PDOP_B=EDOP;     GDOP_B =EDOP;
nGsat_Mtrix = EDOP ;              nBsat_Mtrix = EDOP ;                              nsat_Mtrix = EDOP ;  
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取raw data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
chg = 1 ;                                               % record the change of satellite amounts
Rrg = zeros( 32,2 ) ;                                  Rrb = zeros( 15,2 ) ;
Carrier_G = zeros( 32,2 ) ;                       Carrier_B = zeros( 15,2 ) ;
count_G = zeros( 32,1 ) ;                         count_B = zeros( 15,1 ) ;
Delta_CP_G = zeros( runtime , 32 ) ;      Delta_CP_B = zeros( runtime , 15 ) ;
S_G = zeros( 32,1 ) ;                                  S_B = zeros( 15,1 ) ;
SV_Center = zeros( runtime,3 ) ;
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
%-------------------------------------------------------------------------
mean_delta = zeros( 32,2 ) ;
RMSE = zeros( 32,1 ) ;
ht = zeros( 32,1 ) ;
checkmatrix = zeros( 32,3 ) ;
checkmatrix_B = zeros( 15,3 ) ;
sum0 = 0;
run_trueRange=zeros(32,runtime);
run_ADR=zeros(32,runtime);
run_recbias=zeros(runtime,1);
run_trueRange_B=zeros(15,runtime);
run_ADR_B=zeros(15,runtime);
run_recbias_B=zeros(runtime,1);
cycle_slip_time = zeros(32,runtime);
%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week+10800 ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;
[ constellations_BEIDOU ] = goGNSS.initConstellation( 0 , 0 , 0 , 1 , 0 , 0 ) ;
%   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag,BDS_flag, QZS_flag, SBS_flag);

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
%llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
llh_0 = [ 22.285340 114.161677 3 ] ;       %香港摩天輪
%llh_0 = [-6.367849, 106.833466 3 ] ;        %cibg清真寺
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
%rec_pos_0 = [-1837002.6304  6065627.3599  -716183.2716];
rec_bias_G = 0 ;
rec_bias_B = 0 ;
rec_bias = 0 ;