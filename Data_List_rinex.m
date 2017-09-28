%   此m檔用來讀取定位資料


Data = 0710 ;
no_solid_tide = Data ;
obs_flag = 3 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 60 ;
innovation_window = 80 ;
chenckcount = 0 ;
scale_detect_cycle_slip = 6 ;
initial_time = converge_time + 10 ;
smaller_than_elevation = 0 ;                   %仰角小於elevation即濾掉
interpolation_order = 13 ;                      %內插階數
cut_time = initial_time ;

switch Data
    case{ 0512 }
        addpath 0512
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1948 ;         month = 5 ;      day = 12 ;      day_of_year = 132 ;
        obs_filename = 'rtcm2rnx1.obs' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19485_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg1320.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm1320.17p' ;                                  % 廣播星歷導航資料
    case{ 0710 }
        addpath 20170710
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1957 ;         month = 7 ;      day = 10 ;      day_of_year = 191 ;
        obs_filename = '0710_1100_0.5hr.obs' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19571_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg1910.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm1910.17p' ;                                  % 廣播星歷導航資料
    case{ 0309 }
        addpath 0309
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1939  ;         month = 3 ;      day = 9 ;      day_of_year = 68 ;
        obs_filename = 'rtcm2rnx.obs' ;         % HED
        SP3_filename = 'igu19394_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0680.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0680.17p' ;                                  % 廣播星歷導航資料
    case{ 0623 }
        addpath 0623
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1954 ;         month = 6 ;      day = 23 ;      day_of_year = 174 ;
        obs_filename = '0623_0.5hr.obs' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19545_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg1740.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm1740.17p' ;                                  % 廣播星歷導航資料
    case{ 0829 }
        addpath Data_0829
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 19121 ;         month = 8 ;      day = 29 ;      day_of_year = 242 ;
        obs_filename = 'hkss242hi.16o' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19121_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg2420.16i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm2420.16p' ;                                  % 廣播星歷導航資料
        
    case{ 1024 }
        addpath DG14_1024
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1920 ;         month = 10 ;      day = 24 ;      day_of_year = 298 ;
        obs_filename = '1024.16o' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19201_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg2980.16i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm2980.16p' ;                                  % 廣播星歷導航資料
    case{ 10241 }
        addpath DG14_1024
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1920 ;         month = 10 ;      day = 24 ;      day_of_year = 298 ;
        obs_filename = '1024_0d5ms.16o' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19201_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg2980.16i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm2980.16p' ;                                  % 廣播星歷導航資料
    case{ 0222 }
        addpath 0222
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1937 ;         month = 2 ;      day = 22 ;      day_of_year = 53 ;
        obs_filename = 'pppcmp20170222.17O' ;         % NovAtel
        SP3_filename = 'igu19373_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0530.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'pppcmp20170222.17N' ;                                  % 廣播星歷導航資料
    case{ 02221 }
        addpath 0222
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1937 ;         month = 2 ;      day = 22 ;      day_of_year = 53 ;
        obs_filename = '0722.17o' ;         % HED GPS/BDS
        SP3_filename = 'igu19373_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0530.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'pppcmp20170222.17N' ;                                  % 廣播星歷導航資料
    case{ 0224 }
        addpath 0224
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1937 ;         month = 2 ;      day = 24 ;      day_of_year = 55 ;
        obs_filename = 'COM4_2017-02-24_16.54.24.obs' ;         % BDS與GPS的數據檔案
        SP3_filename = 'igu19375_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0550.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0550.17p' ;                                  % 廣播星歷導航資料
    case{ 0330 }
        addpath 0330
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1942  ;         month = 3 ;      day = 30 ;      day_of_year = 89 ;
        obs_filename = 'rtcm2rnx.obs' ;         % HED
        SP3_filename = 'igu19424_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg0890.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm0890.17p' ;                                  % 廣播星歷導航資料
    case{ 041020 }
        addpath 0410
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1944  ;         month = 4 ;      day = 10 ;      day_of_year = 100 ;
        obs_filename = '041020.obs' ;         % HED
        SP3_filename = 'igu19441_00.sp3' ;                               % IGS精密星歷資料
        GIM_filename = 'c2pg1000.17i' ;                                   % GIM全球電離層 ionex檔名
        nav_filename = 'brdm1000.17p' ;                                  % 廣播星歷導航資料
end

%-----------讀取接收機資訊------------------------
%[XYZ_station,obs,observablesHeader,measurementsInterval]=readRinex302(obs_filename);
[obs]=readRinexBSvHED(obs_filename);
%[XYZ_station,obs,observablesHeader,measurementsInterval]=readRinexBS_V2(obs_filename);

week_num = obs( : , 1 ) ;
time = obs( : , 2) ;

flag = obs( : , 3) ;
prn = obs(: , 4) ;
pr = obs(: , 5);
ADR = obs( : , 6) ;
pr_rate = obs(: , 7) ;
cnr0 = obs(: , 8) ;

%runtime = round( max( time ) ) - round ( min( time ) ) + 1 ;
runtime = 1000 ;
%%

c = 299792458 ;                                 % 光速(m/s)
L1_freq = 1575.42*10^6 ;                % L1 band frequency(Hz)
lambda_G = c /(1575.42*10^6) ;
lambda_B = c /(1561.098*10^6) ;

