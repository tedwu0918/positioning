tic
format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % 光速(m/s)
lambda = c /(1575.42*10^6) ;

Data = 0629 ;
no_solid_tide = Data ;
obs_flag = 1 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 60 ;
chenckcount = 0 ;
scale_detect_cycle_slip = 3 ;
initial_time = converge_time + 10 ;
%smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉
interpolation_order = 13 ;                      %內插階數
cut_time = initial_time ;
switch Data
    case{ 0629 }
        addpath 0629simulator
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2013 ;     week_num = 1771 ;         month = 12 ;      day = 16 ;      day_of_year = 350 ;
        %obs_filename = 'COM5_2017-06-29_16.04.50.obs' ;         % BDS與GPS的數據檔案
        obs_filename = '1.obs' ;         % BDS與GPS的數據檔案
        nav_filename = 'COM5_2017-06-29_16.04.50.nav' ;                                  % 廣播星歷導航資料
end

%-----------讀取接收機資訊------------------------
%[XYZ_station,obs,observablesHeader,measurementsInterval]=readRinexBS_V2(obs_filename);
[obs]=readRinexBSvHED(obs_filename);

week_num = obs( : , 1 ) ;
time = obs( : , 2) ;
flag = obs( : , 3) ;
prn = obs(: , 4) ;
pr = obs(: , 5);
ADR = obs( : , 6) ;
pr_rate = obs(: , 7) ;
cnr0 = obs(: , 8) ;

runtime = round( max( time ) ) - round ( min( time ) ) + 1 ;
%runtime = 700 ;
%cycle_slip_time = zeros(32,runtime);
ref_pos_0=[ -2144855.42      4397605.31  4078049.85 ];                                %
first_P = 0 ;       %   用來記錄第一次成功定位的時刻數
%-----------initialize--------------------------
GMax = 32 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
EDOP_G=EDOP;                        NDOP_G=EDOP;     VDOP_G=EDOP;     HDOP_G=EDOP;     PDOP_G=EDOP;     GDOP_G=EDOP;
nGsat_Mtrix = EDOP ;             
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
chg = 1 ;                                               % record the change of satellite amounts
Rrg = zeros( 32,2 ) ;                                  
Carrier_G = zeros( 32,2 ) ;                       
count_G = zeros( 32,1 ) ;                         
Delta_CP_G = zeros( runtime , 32 ) ;      
S_G = zeros( 32,1 ) ;                               
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
%-------------------------------------------------------------------------
mean_delta = zeros( 32,2 ) ;
RMSE = zeros( 32,1 ) ;
sum0 = 0;
run_trueRange=zeros(32,runtime);
run_ADR=zeros(32,runtime);
run_recbias=zeros(runtime,1);
%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week+10800 ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC時間(0~86400秒)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % 台灣時區為UTC+8

%------------選擇使用哪種衛星系統-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;
%   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag,BDS_flag, QZS_flag, SBS_flag);

%--------------------------------讀取廣播星歷資料----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage
Data
Rhoc_error = [];
EL_error = [];
CNR_error = [];
%----------------猜測初始位置--------------------------
%llh_0 = [ 25.046751 121.517285 3 ] ;       %台北車站
llh_0 = [ 22.285340 114.161677 3 ] ;       %香港摩天輪
%llh_0 = [-6.367849, 106.833466 3 ] ;        %cibg清真寺
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
%rec_pos_0 = [-1837002.6304  6065627.3599  -716183.2716];
rec_bias_G = 0 ;
rec_bias = 0 ;

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
    ACNR_G = [] ;
    
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
            if prn( count ) < 2000 && prn( count ) ~= 1193 

                prn( count ) = prn( count ) -1000 ;
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count )*lambda ;
                
                %if rem(ADR( count ),0.001) ==0
                    Carrier_G( prn(count) , 2 ) = ADR( count )*lambda;
                    if  obs_flag == 3 && i <= initial_time || obs_flag <= 2
                        GP_id( j , 1) = prn( count ) ;
                        Rho0_Gpr( j , 1 ) = pr( count ) ;
                        CNR_G( j , 1) = cnr0( count ) ;
                        j = j+1 ;
                    elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                        GP_id( j , 1) = prn( count ) ;
                        CNR_G( j , 1) = cnr0( count ) ;
                        
                        j = j+1 ;
                    end
               
               % end
                
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
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code_Rinex( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
        elseif obs_flag == 3
            %[ Carrier_G  , RMSE , mean_delta , ht ] = DCDRM( Carrier_G , Rrg , lambda ,  RMSE , mean_delta , ht , scale_detect_cycle_slip ) ;  %check cycle slip
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
        end
        
        Rr_G = Rho0_G ;
        
        Tg = zeros( nGsat , 1 ) ;
        ig = zeros( nGsat , 1 ) ;
        
        %--------------------得到transmission time ( time_tx )----------------------
        [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
            satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
        

        sat_G = XS_tx ;
        SCBg = zeros( nGsat , 1 ) ;

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
            ACNR_G( l , 1) = CNR_G(a , 1) ;
            Atime_tx( l , 1) = time_tx( a , 1) ;
            runtime_GEL(AGP_id(l,1),i)=AEL_G( l , 1);
            l = l + 1 ;
            %else
            %stEL_G(s,1)= GP_id( a , 1 ) ;
            %s = s + 1 ;
            %end
        end
        
        AnGsat = length( AGP_id ) ;
        Ansat = AnGsat  ;
        
        %-------------------計算covariance----------------------------
        EL =  AEL_G  ;
        %CNR = zeros(Ansat,1) ;
        CNR = [ CNR_G  ] ;
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

        sat_G = XS_tx ;
        SCBg = zeros( nGsat , 1 ) ;
                %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------

                Rhoc_G = ARho0_G - ig - Tg + SCBg  ;
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
                    break
                end
                                %}
            end
            
            if Ansat < 4
                Xxyz( i , : ) = rec_pos_0 ;
            else
                Xxyz( i , : ) = rec_pos ;
            end
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置
%}
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
        %{
        %------------紀錄carrier phase&true range-----------------
        for r=1:AnGsat
            run_ADR(AGP_id(r),i) = Carrier_G( AGP_id(r) , 2 ) + ig(r , 1) - Tg(r , 1) + SCBg(r , 1) + relative_IGS(r , 1) - rec_bias ;
            run_trueRange(AGP_id(r),i) = ARr_G(r) ;
        end
        run_recbias(i)= rec_bias;
        %}
        %------------紀錄carrier phase&true range-----------------
        if i > initial_time
            [ sss , ARr_G ] = Fix_Earth_Rotation( sat_G , ref_pos_0 ) ;
            Rhoc_error = [ Rhoc_error ; (Rhoc - ARr_G-rec_bias)];
            EL_error =[EL_error ; EL] ;
            CNR_error = [CNR_error ;CNR] ;
        end
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
if obs_flag ==3
    CutTimeArray = 1 : cut_time ;
    Xxyz(CutTimeArray , :) = [] ;
    nsat_Mtrix(CutTimeArray) = [] ;
    PDOP(CutTimeArray) = [] ;
    sampling_time_end = runtime - cut_time ;
    TimeArray = first_P : sampling_time_end ;
    %effTimeArray_error = round(length(TimeArray)/100) -1 + first_P : sampling_time_end ;
    effTimeArray = TimeArray;
elseif obs_flag == 1
    sampling_time_end = runtime ;
    effTimeArray = 1 : sampling_time_end ;
end

figure( 11 ) ;
scatter(EL_error,CNR_error );
grid on ;
xlabel( 'EL' ) ;
ylabel( 'CNR' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

CNR_error = round(CNR_error) ;
[CNR_error , ic] = sort(CNR_error) ;
CRhoc_error = Rhoc_error(ic) ;

EL_error = round(EL_error) ;
[EL_error , id] = sort(EL_error) ;
ERhoc_error = Rhoc_error(id) ;

cref = CNR_error(1) ;
Eref = EL_error(1) ;

ccalculate = [] ;
Ecalculate = [] ;
n=1;x_matrix = [] ;mean_matrix = [] ;std_matrix = [];
ss = 1 ;Ex_matrix = [] ;Emean_matrix = [] ;Estd_matrix = [];
for i = 1 : length(CNR_error)
    if CNR_error(i)~= cref
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        cref =  CNR_error(i);
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
    end
      if EL_error(i)~= Eref
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Eref =  EL_error(i);
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
    ccalculate(n ,1) = CRhoc_error(i) ;
    Ecalculate(ss,1) = ERhoc_error(i) ;
    n = n+1;
    ss = ss + 1 ;
    
    
        if i==length(CNR_error)
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        cref =  CNR_error(i);
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
    end
      if i==length(CNR_error)
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Eref =  EL_error(i);
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
end
%}

figure( 11 ) ;
scatter(CNR_error , CRhoc_error );
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

figure( 12 ) ;
scatter(EL_error , ERhoc_error );
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('EL-CSC error');
print( '-dpng',  'EL-CSC error' , '-r600' ) ;


figure( 101 ) ;
plot(x_matrix,std_matrix);
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'variance(m)' ) ;
grid on;                        hold on;
title('CNR variance');
print( '-dpng',  'CNR variance' , '-r600' ) ;


figure( 102 ) ;
plot(Ex_matrix,Estd_matrix);
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'variance(m)' ) ;
grid on;                        hold on;
title('EL variance');
print( '-dpng',  'EL variance' , '-r600' ) ;


%mean Plot

 figure( 103 ) ;
bar(x_matrix,mean_matrix);
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('CNR mean error');
print( '-dpng',  'CNR mean error' , '-r600' ) ;


figure( 104 ) ;
bar(Ex_matrix,Emean_matrix);
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('EL mean error');
print( '-dpng',  'EL mean error' , '-r600' ) ;


[ Data_start_end_obs ] = { Data sampling_time_end obs  }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 )       %計算STD
%rec_pos_act = lla2ecef([22.4259617417365 , 114.211355684335 ,53.2231195336208], 'WGS84' )   %華大測站參考位置
Xenu = zeros( runtime , 3 ) ;
for i = 1 : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    sum0 = sum0 + Xenu( i , : ).^2 ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( sum0/runtime ) ;
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;

%mean_errorE = mean(Xenu(effTimeArray,1))
%mean_errorN = mean(Xenu(effTimeArray,2))
%mean_errorU = mean(Xenu(effTimeArray,3))
%mean_errorH = sqrt(mean_errorE^2 + mean_errorN^2)
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
%axis( [-1.5 , 1.5 ,-1.5 , 1.5 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
%axis( [xlim,-2,1] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-2,2] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-5,5] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI

%%

%-------------------------------MSE----------------------------------------------------------
sum0 = 0;
rec_pos_act = [-2144855.42      4397605.31  4078049.85];

Xenu = zeros( runtime , 3 ) ;
for i = 1 : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    sum0 = sum0 + Xenu( i , : ).^2 ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( sum0/runtime ) ;
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;

%mean_errorE = mean(Xenu(effTimeArray,1))
%mean_errorN = mean(Xenu(effTimeArray,2))
%mean_errorU = mean(Xenu(effTimeArray,3))
%mean_errorH = sqrt(mean_errorE^2 + mean_errorN^2)
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
STD = [ std_enu' ; std_en' ]


toc