%   此m檔用來讀取定位資料
%   所需用到的變數為Data

Data = 0512  ;
no_solid_tide = Data ;
obs_flag = 1 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 900 ;
initial_time = converge_time + 20 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
G_B_condition = 1.05 ;
interpolation_order = 13 ;                      %內插階數
Remark = '' ;

switch Data
    case { 0322 }
        addpath Data_ntou0322 ;
        %rec_pos_act = lla2ecef( [ 25.019672 121.540897 58.5 ] , 'WGS84' ) ;
        year = 2016 ;      month = 3 ;     day = 22 ;       day_of_year = 82 ;
        sampling_time_start = 20 ;                       %sampling_time_end = 10991 ;                  solid_filename = 'solid_1117.txt' ;
        obs_filename = 'ntou0322_16000.csv' ;                          SP3_filename = 'igu18892_00.sp3' ;
        GIM_filename = 'c2pg0820.16i' ;                     nav_filename = 'brdm0820.16p' ;
    case { 0323 }
        addpath Data_ntou0323 ;
        %rec_pos_act = lla2ecef( [ 25.019672 121.540897 58.5 ] , 'WGS84' ) ;
        year = 2016 ;      month = 3 ;     day = 23 ;       day_of_year = 83 ;
        sampling_time_start = 20 ;                       %sampling_time_end = 10991 ;                  solid_filename = 'solid_1117.txt' ;
        obs_filename = 'ntou0323_10000.csv' ;                          SP3_filename = 'igu18893_00.sp3' ;
        GIM_filename = 'c2pg0830.16i' ;                     nav_filename = 'brdm0830.16p' ;
    case { 0512 }
        addpath Data_0512 ;
        year = 2016 ;      month = 5 ;     day = 12 ;       day_of_year = 133 ;
        sampling_time_start = 40 ;                       %sampling_time_end = 10991 ;                  solid_filename = 'solid_1117.txt' ;
        obs_filename = 'ntou0512.csv' ;                          SP3_filename = 'igu18964_00.sp3' ;
        GIM_filename = 'c2pg1330.16i' ;                     nav_filename = 'brdm1330.16p' ;
    case { 0516 }
        addpath Data_0516 ;
        year = 2016 ;      month = 5 ;     day = 16 ;       day_of_year = 137 ;
        sampling_time_start = 20 ;                       %sampling_time_end = 10991 ;                  solid_filename = 'solid_1117.txt' ;
        obs_filename = '0516.csv' ;                          SP3_filename = 'igu18971_00.sp3' ;
        GIM_filename = 'c2pg1370.16i' ;                     nav_filename = 'brdm1370.16p' ;
    case { 0530 }
        addpath Data_0530 ;
        year = 2016 ;      month = 5 ;     day = 30 ;       day_of_year = 151 ;         sampling_time_start = 20 ;
        obs_filename = 'HED_0530.csv' ;                     SP3_filename = 'igu18991_00.sp3' ;
        GIM_filename = 'c2pg1510.16i' ;                     nav_filename = 'brdm1510.16p' ;
    case { 06291 }
        addpath Data_0629 ;
        year = 2016 ;      month = 6 ;     day = 29 ;       day_of_year = 181 ;         sampling_time_start = 20 ;
        obs_filename = '0629_1044-1047.csv' ;                     SP3_filename = 'igu19033_00.sp3' ;
        GIM_filename = 'c2pg1810.16i' ;                     nav_filename = 'brdm1810.16p' ;
    case { 06292 }
        addpath Data_0629 ;
        year = 2016 ;      month = 6 ;     day = 29 ;       day_of_year = 181 ;         sampling_time_start = 20 ;
        obs_filename = '0629_1049-1052.csv' ;                     SP3_filename = 'igu19033_00.sp3' ;
        GIM_filename = 'c2pg1810.16i' ;                     nav_filename = 'brdm1810.16p' ;
    case { 0712 }
        addpath Data_0712 ;
        year = 2016 ;      month = 7 ;     day = 12 ;       day_of_year = 194 ;         sampling_time_start = 20 ;
        obs_filename = '0712.csv' ;                     SP3_filename = 'igu19052_00.sp3' ;
        GIM_filename = 'c2pg1940.16i' ;                     nav_filename = 'brdm1940.16p' ;
    case { 0824 }
        addpath Data_0824 ;
        year = 2016 ;      month = 8 ;     day = 24 ;       day_of_year = 237 ;         sampling_time_start = 20 ;
        obs_filename = '0824.csv' ;                     SP3_filename = 'igu19113_00.sp3' ;
        GIM_filename = 'c2pg2370.16i' ;                     nav_filename = 'brdm2370.16p' ;
    case { 0920 }
        addpath Data_0920 ;
        year = 2016 ;      month = 9 ;     day = 20 ;       day_of_year = 264 ;         sampling_time_start = 20 ;
        obs_filename = '0920.csv' ;                     SP3_filename = 'igu19152_00.sp3' ;
        GIM_filename = 'c2pg2640.16i' ;                     nav_filename = 'brdm2640.16p' ;
    case { 12291 }
        addpath Data_1229 ;
        year = 2016 ;      month =12 ;     day = 29 ;       day_of_year = 364 ;         sampling_time_start = 20 ;
        obs_filename = '1229_1.csv' ;                     SP3_filename = 'igu19294_00.sp3' ;
        GIM_filename = 'c2pg3640.16i' ;                     nav_filename = 'brdm3640.16p' ;
    case { 12292 }
        addpath Data_1229 ;
        year = 2016 ;      month =12 ;     day = 29 ;       day_of_year = 364 ;         sampling_time_start = 20 ;
        obs_filename = '1229_2.csv' ;                     SP3_filename = 'igu19294_00.sp3' ;
        GIM_filename = 'c2pg3640.16i' ;                     nav_filename = 'brdm3640.16p' ;
    case { 12293 }
        addpath Data_1229 ;
        year = 2016 ;      month =12 ;     day = 29 ;       day_of_year = 364 ;         sampling_time_start = 20 ;
        obs_filename = '1229_3.csv' ;                     SP3_filename = 'igu19294_00.sp3' ;
        GIM_filename = 'c2pg3640.16i' ;                     nav_filename = 'brdm3640.16p' ;
end

c = 299792458 ;                                 % 光速(m/s)
L1_freq = 1575.42*10^6 ;                % L1 band frequency(Hz)
lambda_G = c /(1575.42*10^6) ;
lambda_B = c /(1561.098*10^6) ;

%-----------讀取接收機資訊------------------------
[ week_num , hwtime , useless , sat_sys , prn , cnr0 , TOW0 , Elapse_Epoch , Elapse_Code , pr , pr_rate , ADR , epr , sat_EL , sat_AZ , time ] = ...
    textread( obs_filename , '%f %f %c %c %f %f %f %f %f %f %f %f %f %f %f %f' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

%-----------initialize--------------------------
runtime = round( max( time ) ) - round ( min( time ) ) + 1  ;
