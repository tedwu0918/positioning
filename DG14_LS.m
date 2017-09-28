tic

format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % ���t(m/s)
lambda = c /(1575.42*10^6) ;

Data = 01091
obs_flag = 3 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
scale_detect_cycle_slip = 3 ;
converge_time = 60 ;
initial_time = converge_time + 10 ;
chenckcount = 0 ;
smaller_than_elevation = 10 ;                   %�����p��elevation�Y�o��
interpolation_order = 13 ;                      %��������


switch Data
    case{ 1210 }
        addpath DG14_1210
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2015 ;     week_num = 1874 ;         month = 12 ;      day = 10 ;      day_of_year = 344 ;
        obs_filename = '1210.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu18744_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg3440.15i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm3440.15p' ;                                  % �s���P���ɯ���
    case{ 0511 }
        addpath DG14_0511
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1896 ;         month = 5 ;      day = 11 ;      day_of_year = 132 ;
        obs_filename = '0511.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu18963_00.sp3' ;                               % IGS��K�P�����
        %GIM_filename = 'gpsg1320.16i' ;                                  % Broadcast ionex�ɦW
        GIM_filename = 'c2pg1320.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm1320.16p' ;                                  % �s���P���ɯ���
    case{ 0621 }
        addpath DG14_0621
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1902 ;         month = 6 ;      day = 21 ;      day_of_year = 173 ;
        obs_filename = '0621.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19022_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg1730.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm1730.16p' ;                                  % �s���P���ɯ���
    case{ 0622 }
        addpath DG14_0622
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1902 ;         month = 6 ;      day = 22 ;      day_of_year = 174 ;
        obs_filename = '0622.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19023_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg1740.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm1740.16p' ;                                  % �s���P���ɯ���
    case{ 08011 }
        addpath DG14_0801
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1908 ;         month = 8 ;      day = 1 ;      day_of_year = 214 ;
        obs_filename = '0801_1.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19081_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg2140.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm2140.16p' ;                                  % �s���P���ɯ���
    case{ 08012 }
        addpath DG14_0801
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1908 ;         month = 8 ;      day = 1 ;      day_of_year = 214 ;
        obs_filename = '0801_2.csv' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19081_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg2140.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm2140.16p' ;                                  % �s���P���ɯ���
    case{ 1104 }
        addpath DG14_1104
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1921 ;         month = 11 ;      day = 4 ;      day_of_year = 309 ;
        obs_filename = '1104.txt' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19215_18.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg3090.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm3090.16p' ;                                  % �s���P���ɯ���
    case{ 12261 }
        addpath DG14_1226
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1929 ;         month = 12 ;      day = 26 ;      day_of_year = 361 ;
        obs_filename = '1.txt' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19291_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg3610.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm3610.16p' ;                                  % �s���P���ɯ���
    case{ 12262 }
        addpath DG14_1226
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2016 ;     week_num = 1929 ;         month = 12 ;      day = 26 ;      day_of_year = 361 ;
        obs_filename = '2.txt' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19291_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg3610.16i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm3610.16p' ;                                  % �s���P���ɯ���
    case{ 01091 }
        addpath DG14_0109
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1931 ;         month = 1 ;      day = 9 ;      day_of_year = 9 ;
        obs_filename = '_msvgpsuv_com4_0109.csv' ;         % static
        SP3_filename = 'igu19311_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg0090.17i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm0090.17p' ;                                  % �s���P���ɯ���
    case{ 0214 }
        addpath DG14_0214
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1936 ;         month = 2 ;      day = 14 ;      day_of_year = 45 ;
        %obs_filename = '0214_receiver_point.txt' ;         % BDS�PGPS���ƾ��ɮ�
        obs_filename = '0214-_base.txt';
        SP3_filename = 'igu19362_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg0450.17i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm0450.17p' ;                                  % �s���P���ɯ���
        TDD_filename = 'bout.txt'     ;                                      %tripple difference data
    case{ 0118 }
        addpath DG14_0118
        rec_pos_act = lla2ecef( [ 25.019625 121.540913 52 ] , 'WGS84' ) ;
        year = 2017 ;     week_num = 1932 ;         month = 1 ;      day = 18 ;      day_of_year = 18 ;
        obs_filename = '1-18-m2.txt' ;         % BDS�PGPS���ƾ��ɮ�
        SP3_filename = 'igu19323_00.sp3' ;                               % IGS��K�P�����
        GIM_filename = 'c2pg0180.17i' ;                                   % GIM���y�q���h ionex�ɦW
        nav_filename = 'brdm0180.17p' ;                                  % �s���P���ɯ���
        TDD_filename = 'bout182.txt'     ;                                      %tripple difference data
end

%-----------Ū����������T------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:���L�Ĥ@��


time(1:9) = [] ;
prn(1:9) = [] ;
SV_x(1:9) = [] ;
SV_y(1:9) = [] ;
SV_z(1:9) = [] ;
ADR(1:9) = [] ;
pr(1:9) = [] ;
pr_rate(1:9) = [] ;
SV_x_dot(1:9) = [] ;
SV_y_dot(1:9) = [] ;
SV_z_dot(1:9) = [] ;
%-----------initialize--------------------------
SVMax = 32 ;
runtime = round( max( time ) ) - round ( min( time ) ) + 1 ;
runtime = 700 ;
EDOP=zeros( 1 , runtime );     NDOP=EDOP;          VDOP=EDOP;         HDOP=EDOP;          PDOP=EDOP;          GDOP=EDOP;
nsat_Mtrix = EDOP ;
count_times = ( min( time ) ) ;               % ��ƪ��ɨ�
count = 1 ;                                                % �ĴX����(Ū��row data)
solid_tide_error = zeros( runtime , 3 ) ;
Xxyz = zeros( runtime , 3 ) ;
wl_count = 0 ;     %�p��GPS/BDS�v���|�N����
chg = 1 ; % record the change of satellite amounts
PRr = zeros( 32,2 ) ;
Carrier = zeros( 32,2 ) ;
Smooth_Count = zeros( 32,1 ) ;
Smooth_Code = zeros( 32,1 ) ;
mean_delta = zeros( 32,2 ) ;
RMSE = zeros( 32,1 ) ;
ht = zeros( 32,1 ) ;

runtime_GEL=zeros(32,runtime);
runtime_carrier=zeros(32,runtime);
runtime_trueRange=zeros(32,runtime);
runtime_cbias=zeros(runtime);
runtime_error = zeros(32,runtime) ;
runtime_mean_error = zeros(32,runtime) ;

firstpostime = 0 ;
first_pos = 1;
checkmatrix = zeros( 32,3 ) ;

%-----------------set time--------------------------
week = week_num( 1 ) ;                                       % GPS week
time_week = min( time ) ;                         % time of week(start time)
time_week_last_epoch = time_week ;
time_day_UTC = rem( round( time_week ) , 86400 ) ;                  % UTC�ɶ�(0~86400��)
local_time = rem( round( min( time ) ) + 8*3600 , 86400 ) ;        % �x�W�ɰϬ�UTC+8

%------------��ܨϥέ��ؽìP�t��-----------------
[ constellations_GPS ] = goGNSS.initConstellation( 1 , 0 , 0 , 0 , 0 , 0 ) ;

%--------------------------------Ū���s���P�����----------------------------------
[ Eph_tatol , iono_para ] = RINEX_get_nav( nav_filename , constellations_GPS ) ;   % load navigation massage

%----------------------------�o��SP3�ɮת����---------------------------
[ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;

%---------------�qionex�ɮפ��o��VTEC�Plat,lon,time����T----------
[ VTEC_MAP , time_total , lat_total , lon_total ] = ...    % �`�NVTEC����쬰0.1TECU
    GPS_ParseIONEXTEC( GIM_filename ) ;

%----------------�q����l��m--------------------------
llh_0 = [ 25.046751 121.517285 3 ] ;       %�x�_����
rec_pos_0 = lla2ecef( llh_0 , 'WGS84' ) ;
rec_bias = 0 ;
bar = waitbar(0,'Please wait...') ;

for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %�C900��,���s��slagrange�����`�I
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX��(prn*data_num)���x�},�Y �Ĥ@�C���ìPprn1�b13�Ӯɨ誺x��m ,�ĤG�C���ìPprn2�b13�Ӯɨ誺x��m
    
    rec_bias_satpos = 0 ;
    id = [] ;
    Rho0 = [] ;
    Rho0_pr = [] ;
    CNR = [] ; Ds =[];
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
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    
    while judge == 1
        
        if time( count ) == count_times                                                       % all( )==1 , �YAX�x�}�������������s
            if prn( count ) ~= 193  && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999%...
                    %&&(  prn( count ) == 11   || prn( count ) == 9 || prn( count ) == 23 || prn( count ) == 8||prn( count ) == 27  || prn( count ) == 26)
                
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
                    Ds( j , 1 ) = -pr_rate( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && Smooth_Count( prn(count),1 ) > converge_time && Carrier(prn(count),2) ~= 0 && Carrier(prn(count),1) ~= 0
                    id( j , 1) = prn( count ) ;
                    [ EL( j,1) , AZI( j,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , SVxyz ) ;
                    checkmatrix( id( j , 1),1 ) = 1 ;
                    checkmatrix( id( j , 1),3 ) = checkmatrix( id( j , 1),3 ) + 1 ;
                    Ds( j , 1 ) = -pr_rate( count ) ;
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
    
    %----------------�C��(�C�@��)�w��ìP������---------------
    nsat = length( id ) ;
    nsat_Mtrix( i ) = nsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nsat ~= 0
        
        %-----------�o�챵���ɶ�---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %-------------------�p��covariance----------------------------
        CNR = zeros(nsat,1) ;
        [ sat_var_elevation , sat_var_CNR , sat_var_SIGMA , sat_var_CandE ] = weightingfunc( EL , CNR ) ;
        %-------------------����v���x�}------------------------------
        W = diag( ( 1./( sat_var_elevation ) ).^2 ) ;
        %W = eye( nsat ) ;
        if obs_flag == 3 && chenckcount ~=0 && i > initial_time
            [W , checkmatrix ] = adjweightingfunc( W , checkmatrix , id , chenckcount ) ;
        end
        %----------------����[���q------------------
        if obs_flag == 1
            Rho0 = Rho0_pr ;
        elseif obs_flag == 2
            [ Rho0 , Smooth_Code , Smooth_Count ] = Doppler_Smoothed_Code_DG14( PRr , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        elseif obs_flag == 3
          %  [ Carrier , RMSE , mean_delta , ht ] = DCDRM_DG14( Carrier , PRr , lambda ,  RMSE , mean_delta , ht , scale_detect_cycle_slip ) ;  %check cycle slip
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
            %--------------------�o��transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
                satellite_positions( r_gpst , Rr , id , Eph_tatol , [] , [] , Tropo , Iono , rec_bias_satpos ) ;    %�`�N,���N�ɭn�a�J�a�y�ץ��᪺�����Z��(Rr)
            
            %------------�Q�� UNB3 �ҫ��p���y�h�~�t-------------------------
            for n = 1 : nsat
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL(n)*pi/180 ) ;
                Tropo( n , 1 ) = RTROP ;
            end
            
            %------------�ϥ�GIM�p��q���h�~�t--------------------------
            for n = 1 : nsat
                [ iono , VTEC ]=...
                    GIM_corr( AZI(n) , EL(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                Iono( n , 1 ) = iono ;
            end
            
            %------------�ϥ� time_tx �p���K�P����m-----------------
            %--------------initialize----------------
                            sat = zeros( nsat , 3 ) ;
                SCBg = zeros( nsat , 1 ) ;
            %sat = XS_tx ;
            %SCBg = c.*dtS ;
            
            %--------------------�p���K�P���ìP��m�����t---------------------------
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
            
            %------------------�p��GPS�۹�׻~�t-------------------------------
            relative_IGS = zeros( nsat , 1 ) ;
            tgd = 0 ;
            
            for n = 1 : nsat
                icol = find_eph( Eph_tatol , id( n ) , r_gpst ) ;
                Eph = Eph_tatol( : , icol ) ;
                
                [ satp , satv ] = satellite_orbits( time_tx( n ) , Eph , icol , [] ) ;
                
                Sp3_corr = -2*dot( satp , satv )/c ;
                relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
            end
            %}
            %--------------------�ץ������Z��(�U�����Ҭ�����)--------------------------
            % �_�檺SCBb�]�tsat_bias,relative�Pgroup_delay
            Rhoc = Rho0 - Iono - Tropo + SCBg + relative_IGS ;
            
            %---------------------------�ץ��a�y����--------------------------------
            [ Xs , Rr ] = Fix_Earth_Rotation( sat , rec_pos_0 ) ;
            
            %-----------------------�̤p����k�w��-------------------------
            if nsat > 3     %   ��t�Φܤֻݭn4���ìP�өw��
               % [ Delta_Mtrix , Delta_Rho ] = LeastSquare( Xs , rec_pos_0 , Rhoc , W , nsat , rec_bias ) ;
                [ Delta_Mtrix , Delta_Rho, Delta_Ds ] = LeastSquare_PD( Xs , VS_tx  , rec_pos_0 , Rhoc,Ds , W , nsat , rec_bias) ;             
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
            
            if first_pos ==1
                firstpostime = i
                first_pos = 0 ;
            end
        end
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
        rec_pos_0 = Xxyz( i , : ) ;            %   �N�ұo����m��s���U�@�ɨ褧��l��m
        
    end
    
    %--------------����CP&trueRANGE------------
    for r = 1:nsat
 
        runtime_carrier(id(r),i) = Carrier( id(r) , 2 ) + Iono(r) - Tropo(r) + SCBg(r) + relative_IGS (r) - rec_bias;
        runtime_trueRange(id(r),i)=Rr(r);
        
    end
           runtime_cbias(i) = rec_bias;
    
    %----------------�ɶ�+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    if i > 1 && nsat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    
    nsat_prior = nsat ; %   �����W�@��ìP����
    
end
runtime_error = runtime_carrier - runtime_trueRange;
%runtime_mean_error = bsxfun(@minus, runtime_error. - mean(runtime_error,1));
meanerror=mean(runtime_error,2);
runtime_error = bsxfun(@minus , runtime_error , meanerror ) ;
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
firstpostime = 1;
%---------------------��ƨ��ˮɶ�----------------------
if obs_flag == 3  %g�ϥ�CSC�ɪ��e�ɶ��e���C�J�p��
    
    Xxyz(1:initial_time , :) = [] ;
    nsat_Mtrix(1:initial_time ) = [] ;
    PDOP(1:initial_time ) = [] ;
    sampling_time_end = runtime - initial_time ;
    TimeArray = firstpostime : sampling_time_end ;
    effTimeArray_error = TimeArray ;
else
    %Xxyz(1 : 2 , :) = [] ;
    %nsat_Mtrix(1 : 2) = [] ;
    sampling_time_end = runtime  ;
    TimeArray = firstpostime : sampling_time_end ;
    effTimeArray_error = round(length(TimeArray)/100) -1 + firstpostime : sampling_time_end ;
end

Time = length(TimeArray) ;


%---------------------��ܰѼƳ]�w---------------------
[ Data_start_end_obs ] = { Data sampling_time_end obs  }

%--------------------�Nxyz�নenu------------------------
rec_pos_act = mean( Xxyz(TimeArray,:),1 ) ;      %�p��STD
%rec_pos_act = lla2ecef( [25.149442 , 121.777873 , 73.6] , 'WGS84' );
Xenu = zeros( runtime , 3 ) ;
MSE_enu = zeros(1,3);

for i = firstpostime : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    
end
Xenu = Xenu - solid_tide_error ;
for i = firstpostime : sampling_time_end
    MSE_enu = MSE_enu + Xenu( i , : ).^2;
end
%--------------------�p��STD & RMS------------------------
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;
std = sqrt(length(TimeArray) \ MSE_enu) ;
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------�|�ˤ��J(round)����p���I��X��---------------------
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
plot( TimeArray , Xenu( TimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU��m�P�ɶ����Y (GPS)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(TimeArray, Xenu( TimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-4,4] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(TimeArray, Xenu( TimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-20,20] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI


figure( 10 ) ;
plot( TimeArray , runtime_cbias( TimeArray ) , 'b-' ) ;
axis( [xlim,-20,5] ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'receiver bias' )
print( '-dpng',  'receiver bias' , '-r600'  ) ;

toc