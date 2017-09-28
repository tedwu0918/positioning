format long g
clc ;
clear ;
close all ;
addpath Functions

c = 299792458 ;                                 % ���t(m/s)
lambda = c /(1575.42*10^6) ;

Data = 08011
obs_flag = 1 ;                                             % 0=epr, 1=pr, 2=dsc, 3=csc
converge_time = 900 ;
initial_time = converge_time + 20 ;
smaller_than_elevation = 10 ;                   %�����p��elevation�Y�o��
interpolation_order = 13 ;                      %��������

%---------------------TaTi��l�]�w--------------------------------
epk = 1 ;       %   �ֿndata�����
ac_lim = 1       %   �ֿndata����ƪ��W�� ( �T�w��Ƭ�(ac_lim-1),�B�ʬ�Ƭ�1 )
interval = 1        %   �ֿndata�������j�ɨ��
first_P = 0 ;       %   �ΨӰO���Ĥ@�����\�w�쪺�ɨ��
CNR_ac = [] ;
Rho0_ac = [] ;
EL_ac = [] ;
AZI_ac = [] ;
sat_set = [] ;
SCBg_set = [] ;
relative_IGS_set = [] ;

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
end

%-----------Ū����������T------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:���L�Ĥ@��

%-----------initialize--------------------------
SVMax = 32 ;
runtime = round( max( time ) ) - round ( min( time ) ) ;
%runtime = 3600 ;
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
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    
    while judge == 1
        
        if time( count ) == count_times                                                       % all( )==1 , �YAX�x�}�������������s
            if prn( count ) ~= 193 && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999
                
                SVxyz(1,1) = SV_x(count) ;
                SVxyz(1,2) = SV_y(count) ;
                SVxyz(1,3) = SV_z(count) ;
                 
                PR( prn(count) , 1 ) = pr( count ) ;
                PRr( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier( prn(count) , 2 ) = ADR( count ) ;          %* lambda ;       
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    id( j , 1) = prn( count ) ;
                    [ EL( j,1) , AZI( j,1) ] = Calc_Azimuth_Elevation( rec_pos_act , SVxyz ) ;
                    Rho0_pr( j , 1 ) = pr( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && Smooth_Count( prn(count),1 ) > converge_time && Carrier(prn(count),2) ~= 0 && Carrier(prn(count),1) ~= 0
                    id( j , 1) = prn( count ) ;      
                    [ EL( j,1) , AZI( j,1) ] = Calc_Azimuth_Elevation( rec_pos_act , SVxyz ) ;
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
                
        %----------------����[���q------------------
        if obs_flag == 1
            Rho0 = Rho0_pr ;
        elseif obs_flag == 2
            [ Rho0 , Smooth_Code , Smooth_Count ] = Doppler_Smoothed_Code( PRr , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        elseif obs_flag == 3
            [ Rho0 , Smooth_Code , Smooth_Count ] = Carrier_Smoothed_Code( Carrier , PR , Smooth_Code , id , Smooth_Count , SVMax ) ;
        end
        
        %-------------------�ǤJcurrent data----------------------------
        nID( epk ) = nsat ;
        CNR_ac = [ CNR_ac ; CNR ] ;
        Rho0_ac = [ Rho0_ac ; Rho0 ] ;
        EL_ac = [ EL_ac ; EL ] ;
        AZI_ac = [ AZI_ac ; AZI ] ;
        
        %-------------------�p��covariance----------------------------
        [ sat_var_elevation , sat_var_CNR , sat_var_SIGMA , sat_var_CandE ] = weightingfunc( EL_ac , AZI_ac ) ;
        %-------------------����v���x�}------------------------------
        W = diag( ( 1./( sat_var_elevation ) ).^2 ) ;
        %W = eye( nsat ) ;
        
        Rr = Rho0 ;
        clear Delta_Mtrix 
        Delta_Mtrix(1:4) = 1 ;
        
        while norm( Delta_Mtrix(end-3:end-1) ) > 0.01
            
            llh_0 = xyz2llh( rec_pos_0 ) ;
            lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
            lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
            height_rec = llh_0( 3 ) ;                  % unit : meter         
            
            Tropo = zeros( sum(nID) , 1 ) ;
            Iono = zeros( sum(nID) , 1 ) ;
            %--------------------�o��transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
                satellite_positions( r_gpst , Rr , id , Eph_tatol , [] , [] , Tropo , Iono , rec_bias_satpos ) ;    %�`�N,���N�ɭn�a�J�a�y�ץ��᪺�����Z��(Rr)
            
            %------------�Q�� UNB3 �ҫ��p���y�h�~�t-------------------------
            for n = 1 : sum(nID)
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec*pi/180 , height_rec , day_of_year , EL_ac(n)*pi/180 ) ;
                Tropo( n , 1 ) = RTROP ;
            end
            
            %------------�ϥ�GIM�p��q���h�~�t--------------------------
            for n = 1 : sum(nID)
                [ iono , VTEC ]=...
                    GIM_corr( AZI_ac(n) , EL_ac(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                Iono( n , 1 ) = iono ;
            end
            
            %------------�ϥ� time_tx �p���K�P����m-----------------
            %--------------initialize----------------
            sat = zeros( nsat , 3 ) ;
            SCBg = zeros( nsat , 1 ) ;
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
            %----------------------�ϥ�saved data + current data---------------------------
            sat_ac = [ sat_set ; sat ] ;
            
            %---------------------------�ץ��a�y����--------------------------------
            [ Xs , Rr_temp ] = Fix_Earth_Rotation( sat_ac , rec_pos_0 ) ;
            Rr = Rr_temp( end-nsat+1:end ) ;
            
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
            %----------------------�ϥ�saved data + current data---------------------------
            SCBgac = [ SCBg_set ; SCBg ] ;
            relative_IGSac = [ relative_IGS_set ; relative_IGS ] ;
            
            %--------------------�ץ������Z��(�U�����Ҭ�����)--------------------------
            % �_�檺SCBb�]�tsat_bias,relative�Pgroup_delay
            Rhoc = Rho0_ac - Iono - Tropo + SCBgac + relative_IGSac ;
            
            %-----------------------�̤p����k�w��-------------------------   
            if nsat > 3     %   ��t�Φܤֻݭn4���ìP�өw��
                      
                 [ Delta_Mtrix , Delta_Rho , G ] = LeastSquare_tati_K( Xs , rec_pos_0 , Rhoc , W , nID ) ;
                 
                rec_pos = rec_pos_0 + Delta_Mtrix(end-3:end-1)' ;
                rec_bias = Delta_Mtrix(end-1) ;
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
        rec_pos_0 = Xxyz( i , : ) ;            %   �N�ұo����m��s���U�@�ɨ褧��l��m        
        
        if nsat > 3 && first_P == 0 && all( Xxyz( i,: ) ) == 1  %�Ĥ@���w�즨�\��,������U�ɨ�
            first_P = i ;
        end
         cond = rem( i - first_P , interval ) ;
        %----------------------------�ܮw�x�}(�x�s�®ɨ褧���)---------------------------------------
        if cond ~= 0            % ���ŦX�϶��ɨ�,�M��current data
            Rho0_ac( end-nsat+1:end ) = [] ;
            EL_ac( end-nsat+1:end ) = [] ;
            AZI_ac( end-nsat+1:end ) = [] ;      
        elseif  cond == 0 && epk < ac_lim            % �ŦX�϶��ɨ�B�ֿn�ɨ�q�p��ֿn�W��
            sat_set = sat_ac ;
            SCBg_set = SCBgac ;
            relative_IGS_set = relative_IGSac ;
            epk = epk + 1 ;
            first_P = i ;
        elseif  cond == 0 && epk == ac_lim            % �ŦX�϶��ɨ�B�ֿn�ɨ�q����ֿn�W��,�M��oldest data            
            Rho0_ac( end-nID(epk)+1:end ) = [] ;%end-nID(epk)+1:end
            EL_ac( end-nID(epk)+1:end ) = [] ;
            AZI_ac( end-nID(epk)+1:end ) = [] ;
            sat_set = sat_ac( 1:end-nID(epk),: ) ;
            SCBg_set = SCBgac( 1:end-nID(epk),: ) ;
            relative_IGS_set = relative_IGSac( 1:end-nID(epk),: ) ;
            nID(epk) = [] ;
            first_P = i ;
        end
        
    end
    
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

close( bar ) ;
switch obs_flag
    case { 0 }
        obs = 'EPR' ;
    case { 1 }
        obs = 'PR' ;
    case { 2 }
        obs = 'DSC' ;
    case { 3 }
        obs = 'CSC' ;
end
%---------------------��ƨ��ˮɶ�----------------------
sampling_time_end = runtime ;
effTimeArray = 10 : sampling_time_end ;

%---------------------��ܰѼƳ]�w---------------------
[ Data_start_end_obs_elevation_GBcond ] = { Data sampling_time_end obs  }

%--------------------�Nxyz�নenu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;      %�p��STD

Xenu = zeros( runtime , 3 ) ;
for i = 1 : runtime 
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------�p��STD & RMS------------------------
std = sqrt( var( Xenu(effTimeArray,:))  ) ;
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------�|�ˤ��J(round)����p���I��X��---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ; 
STD = [ std_enu' ; std_en' ] 

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
xlabel( 'East (m)' ) ;
ylabel( 'North (m)') ;
axis equal ;
%axis( [-0.8 , 0.8 , -0.5 , 1.2 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

