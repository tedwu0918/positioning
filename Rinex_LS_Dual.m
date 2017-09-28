tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List_rinex                  %   �qData_List_rinex.m��Ū�����

cycle_slip_time = zeros(32,runtime);
cycle_slip_time_B = zeros(15,runtime);

Preparation               %   �e�m�@�~
Data


%%
%error refence
Rhoc_error_G = [];
EL_error_G = [];
CNR_error_G = [];
Rhoc_error_B = [];
EL_error_B = [];
CNR_error_B = [];
ref_pos_0 = lla2ecef([22.4259617417365 , 114.211355684335 ,53.2231195336208], 'WGS84' ) %ref

%%
first_P = 0 ;       %   �ΨӰO���Ĥ@�����\�w�쪺�ɨ��
bar = waitbar(0,'Please wait...');
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %�C900��,���s��slagrange�����`�I
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX��(prn*data_num)���x�},�Y �Ĥ@�C���ìPprn1�b13�Ӯɨ誺x��m ,�ĤG�C���ìPprn2�b13�Ӯɨ誺x��m
    hel_count=0;
    wl_count = 0 ;                                      %�p��GPS/BDS�v���|�N����
    id_temp_G = [] ;                                       id_temp_B = [] ;
    GP_id = [] ;                                                BD_id = [] ;
    Rho0_G = [] ;                                            Rho0_B = [] ;
    Rho0_Gpr = [] ;                                        Rho0_Bpr = [] ;
    CNR_G = [] ;                                             CNR_B = [] ;
    EL_G = [] ;                                                  EL_B = [] ;
    AZ_G = [] ;                                                 AZ_B = [] ;
    Rg = zeros( 32,1 ) ;                                    Rb = zeros( 15,1 ) ;
    Rrg( :,1 ) = Rrg( :,2 ) ;                                 Rrb( :,1 ) = Rrb( :,2 ) ;
    Rrg( :,2 ) = 0 ;                                             Rrb( :,2 ) = 0 ;
    Carrier_G( :,1 ) = Carrier_G( :,2 ) ;            Carrier_B( :,1 ) = Carrier_B( :,2 ) ;
    Carrier_G( :,2 ) = 0 ;                                  Carrier_B( :,2 ) = 0 ;
    
    AEL_G = [] ;            AAZ_G = [] ;
    ARr_G = [] ;            ARho0_G = [] ;
    ASCBg = [] ;            Asat_G = [] ;
    AGP_id = [] ;            ACNR_G = [] ;
    ACNR_G = [] ;
    Atime_tx = [] ;
    AEL_B = [] ;            AAZ_B = [] ;
    ARr_B = [] ;            ARho0_B = [] ;
    ASCBB = [] ;            Asat_B = [] ;
    ABD_id = [] ;            ACNR_B = [] ;
    Atime_tx_B = [] ;
    ACNR_B = [] ;
    stEL_G = [] ;
    m = 1 ;
    l = 1 ;
    s = 1 ;
    first_in = 0 ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    k = 1 ;                           %�C�@����U��BDS����ƨ̧ǧ�i��
    
    while judge == 1
        
        
        
        if time( count ) == count_times                                                       % all( )==1 , �YAX�x�}�������������s
            if prn( count ) < 2000 && prn( count ) ~= 1193 && all( AX( prn(count)-1000,: ) ) == 1 && AS( prn(count)-1000,1 ) ~= 999999.999999%...
                    %&&( prn( count ) == 1003  || prn( count ) == 1014  || prn( count ) == 1016 || prn( count ) == 1021|| prn( count ) == 1023 || prn( count ) == 1026|| prn( count ) == 1031 || prn( count ) == 1027);
                %&&( prn( count ) == 1003  || prn( count ) == 1014  || prn( count ) == 1022 || prn( count ) == 1026|| prn( count ) == 1031 || prn( count ) == 1032);
                prn( count ) = prn( count ) -1000 ;
                
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count )*lambda_G ;
                if rem(ADR( count ),0.001) ==0
                    Carrier_G( prn(count) , 2 ) = ADR( count )*lambda_G;
                    checkmatrix( prn(count),1 ) = 1 ;
                    checkmatrix( prn(count),3 ) = checkmatrix( prn(count),3 ) + 1 ;
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
                end
            elseif prn( count ) > 5000 && prn( count ) ~= 1193 && all( AX( prn(count)-5000,: ) ) == 1 && AS( prn(count)-5000,1 ) ~= 999999.999999%...
                    % &&( prn( count ) == 5006  || prn( count ) == 5007  || prn( count ) == 5008 || prn( count ) == 5009|| prn( count ) == 5010 || prn( count ) == 5013);
               prn( count ) = prn( count ) -5000 ;
                
                
                Rb( prn(count) , 1 ) = pr( count ) ;
                Rrb( prn(count) , 2 ) = pr_rate( count )*lambda_B ;
                if rem(ADR( count ),0.001) ==0
                    Carrier_B( prn(count) , 2 ) = ADR( count )*lambda_B;
                    checkmatrix_B( prn(count),1 ) = 1 ;
                    checkmatrix_B( prn(count),3 ) = checkmatrix_B( prn(count),3 ) + 1 ;
                    if  obs_flag == 3 && i <= initial_time || obs_flag <= 2
                        BD_id( k , 1) = prn( count ) ;
                        Rho0_Bpr( k , 1 ) = pr( count ) ;
                        CNR_B( k , 1) = cnr0( count ) ;
                        k = k+1 ;
                    elseif obs_flag == 3 && i > initial_time && count_B( prn(count),1 ) > converge_time && Carrier_B(prn(count),2) ~= 0 && Carrier_B(prn(count),1) ~= 0
                        BD_id( k , 1) = prn( count ) ;
                        CNR_B( k , 1) = cnr0( count ) ;                       
                        k = k+1 ;
                    end
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
    nGsat = length( GP_id ) ;
    nBsat = length( BD_id ) ;
    nsat = nGsat + nBsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nGsat ~= 0 && nBsat ~= 0
       
        %-----------�o�챵���ɶ�---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %----------------����[���q------------------
        if obs_flag == 1
            Rho0_G = Rho0_Gpr ;
            Rho0_B = Rho0_Bpr ;
        elseif obs_flag == 2
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
            [ Rho0_B , S_B , count_B ] = Doppler_Smoothed_Code( Rrb , Rb , S_B , BD_id , count_B , BMax ) ;
        elseif obs_flag == 3
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
            [ Rho0_B , S_B , count_B ] = HED_Carrier_Smoothed_Code( Carrier_B , Rb , S_B , BD_id , count_B , BMax , lambda_B ) ;
        end
        Tg = zeros( nGsat , 1 ) ;
        ig = zeros( nGsat , 1 ) ;
        Tb = zeros( nBsat , 1 ) ;
        ib =  zeros( nBsat , 1 ) ;
        Rr_G = Rho0_G ;
        Rr_B = Rho0_B ;
        %--------------------�o��transmission time ( time_tx )----------------------
        [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
            satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
        [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
            satellite_positions( r_gpst , Rr_B , BD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
        
        %------------�ϥ� time_tx �p���K�P����m-----------------
        %--------------initialize----------------
        sat_G = zeros( nGsat , 3 ) ;
        sat_B = XS_tx_B ;
        %--------------------�p���K�P���ìP��m�����t---------------------------
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
        end
        
        %-------------------�p�����----------------------------------
        for a= 1 : nGsat
            [ EL_G( a , 1 ) , AZ_G( a,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , sat_G( a , : ) ) ;
            %if EL_G( a , 1 ) > smaller_than_elevation
            AEL_G( l , 1) = EL_G( a , 1 ) ;
            AAZ_G( l , 1) = AZ_G( a , 1 ) ;
            ARr_G( l , 1) = Rr_G( a , 1 ) ;
            ARho0_G( l , 1) = ARr_G( l , 1) ;
            AGP_id( l , 1) = GP_id(a , 1) ;
            ACNR_G( l , 1) = CNR_G(a , 1) ;
            l = l + 1 ;
            %else
            %stEL_G(s,1)= GP_id( a , 1 ) ;
            %s = s + 1 ;
            %end
        end
        l = 1 ;
        for a= 1 : nBsat
            [ EL_B( a , 1 ) , AZ_B( a,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , sat_B( a , : ) ) ;
            %if EL_G( a , 1 ) > smaller_than_elevation
            AEL_B( l , 1) = EL_B( a , 1 ) ;
            AAZ_B( l , 1) = AZ_B( a , 1 ) ;
            ARr_B( l , 1) = Rr_B( a , 1 ) ;
            ARho0_B( l , 1) = ARr_B( l , 1) ;
            ABD_id( l , 1) = BD_id(a , 1) ;
            ACNR_B( l , 1) = CNR_B(a , 1) ;
            l = l + 1 ;
            %else
            %stEL_G(s,1)= GP_id( a , 1 ) ;
            %s = s + 1 ;
            %end
        end
        AnGsat = length( AGP_id ) ;
        AnBsat = length( ABD_id ) ;
        Ansat = AnGsat + AnBsat ;
        
        %-------------------�p��covariance----------------------------
        EL = [ AEL_G ; AEL_B ] ;
        CNR = [ ACNR_G ; ACNR_B ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        
        %-------------------����v���x�}------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( Ansat ) ;
        %{
        if obs_flag == 3 && chenckcount ~=0 && i > initial_time
            [W , checkmatrix ] = adjweightingfunc( W , checkmatrix , AnGsat , chenckcount ) ;
        end
        %}
        if AnGsat ~= 0 && AnBsat ~= 0
            
            Delta_Mtrix( 1:5 ) = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 0.01
                
                llh_0 = xyz2llh( rec_pos_0 ) ;
                lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
                lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
                height_rec = llh_0( 3 ) ;                  % unit : meter
                
                Tg = zeros( AnGsat , 1 ) ;
                ig = zeros( AnGsat , 1 ) ;
                Tb = zeros( AnBsat , 1 ) ;
                ib =  zeros( AnBsat , 1 ) ;
                
                %--------------------�o��transmission time ( time_tx )----------------------
                [ XS , dtS , XS_tx , VS_tx , Atime_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
                    satellite_positions( r_gpst , ARr_G , AGP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
                [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
                    satellite_positions( r_gpst , ARr_B , ABD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
                
                %------------�Q�� UNB3 �ҫ��p���y�h�~�t-------------------------
                for n = 1 : AnGsat
                    [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , AEL_G(n) * pi/180 ) ;
                    Tg( n , 1 ) = RTROP ;
                end
                for n = 1 : AnBsat
                    [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , AEL_B(n) * pi/180 ) ;
                    Tb( n , 1 ) = RTROP_B ;
                end
                %------------�ϥ�GIM�p��q���h�~�t--------------------------
                for n = 1 : AnGsat
                    [ iono , VTEC ]=...
                        GIM_corr( AAZ_G(n) , AEL_G(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                    ig( n , 1 ) = iono ;
                end
                for n = 1 : nBsat
                    [ iono_B , VTEC_B ]=...
                        GIM_corr( AAZ_B(n) , AEL_B(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                    ib( n , 1 ) = iono_B ;
                end
                %------------�ϥ� time_tx �p���K�P����m-----------------
                %--------------initialize----------------
                sat_G = zeros( AnGsat , 3 ) ;
                SCBg = zeros( AnGsat , 1 ) ;
                
                %--------------------�p���K�P���ìP��m�����t---------------------------
                for k = 1 : AnGsat
                    sat_G( k , 1 ) = LagrangeInter( Data_time , AX( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    sat_G( k , 2 ) = LagrangeInter( Data_time , AY( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    sat_G( k , 3 ) = LagrangeInter( Data_time , AZ( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    for t = 1 : length( AX( AGP_id(k),: ) )
                        if Data_time( t ) - Atime_tx( k ) > 0
                            Interval = t-1 : t ;
                            break
                        end
                    end
                    ephemeris_SCB = AS( AGP_id(k),: ) ; %SCB=satellite clock bias
                    SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , Atime_tx( k ) ) * c * 1e-6 ;
                end
                
                
                %-------------------�ϥ�stallie_position------------------------
                sat_B = XS_tx_B ;
                SCBb = c .* dtS_B ;
                %------------------�p��GPS�۹�׻~�t-------------------------------
                relative_IGS = zeros( AnGsat , 1 ) ;
                tgd = 0 ;
                
                for n = 1 : AnGsat
                    icol = find_eph( Eph_tatol , AGP_id( n ) , r_gpst ) ;
                    Eph = Eph_tatol( : , icol ) ;
                    
                    [ satp , satv ] = satellite_orbits( Atime_tx( n ) , Eph , icol , [] ) ;
                    
                    Sp3_corr = -2*dot( satp , satv )/c ;
                    relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
                end
                
                %--------------------�ץ������Z��(�U�����Ҭ�����)--------------------------
                % �_�檺SCBb�]�tsat_bias,relative�Pgroup_delay
                Rhoc_G = ARho0_G - ig - Tg + SCBg + relative_IGS ;
                Rhoc_B = ARho0_B - ib - Tb + SCBb ;
                Rhoc = [ Rhoc_G ; Rhoc_B ] ;
                
                
                %---------------------------�ץ��a�y����--------------------------------
                [ Xs_G , ARr_G ] = Fix_Earth_Rotation( sat_G , rec_pos_0 ) ;
                [ Xs_B , ARr_B ] = Fix_Earth_Rotation( sat_B , rec_pos_0 ) ;
                Xs = [ Xs_G ; Xs_B ] ;
                %-----------------------�̤p����k�w��-------------------------
                if Ansat > 4
                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                    wl_count = 0;
                    while  ( abs(V_hyber(1)/V_hyber(2)) > G_B_condition || abs(V_hyber(2)/V_hyber(1)) > G_B_condition )
                        W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat , 1:nGsat) , ( 1/V_hyber(2) ) * W(nGsat+1:nsat , nGsat+1:nsat) ) ;
                        
                        [ Delta_Mtrix , Delta_Rho , V_hyber ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                        wl_count = wl_count + 1 ;
                        if wl_count >1000
                            hel_count = hel_count + 1 ;
                            break
                        end
                    end
                    
                    rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias_G = rec_bias_G + Delta_Mtrix(4) ;
                    rec_bias_B = rec_bias_B + Delta_Mtrix(5) ;
                    rec_pos_0 = rec_pos ;
                    
                else
                    break
                end
                if hel_count >= 100
                    i
                    norm( Delta_Mtrix(1:3) )
                    AGP_id
                    ABD_id
                    break
                end
            end
            
            if Ansat < 5
                Xxyz( i , : ) = rec_pos_0 ;
            else
                Xxyz( i , : ) = rec_pos ;
            end
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   �N�ұo����m��s���U�@�ɨ褧��l��m
            %--------------�N�C�����ìP��CSC��ƲM��----------
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
        %------------����carrier phase&true range-----------------
        if i > initial_time
            [ sss , ARr_G ] = Fix_Earth_Rotation( sat_G , ref_pos_0 ) ;
            [ Xs_B , ARr_B ] = Fix_Earth_Rotation( sat_B , ref_pos_0 ) ;
            Rhoc_error_G = [ Rhoc_error_G ; (Rhoc_G - ARr_G-rec_bias_G) ];
             Rhoc_error_B = [ Rhoc_error_B ; (Rhoc_B - ARr_B-rec_bias_B) ];           
            EL_error_G =[EL_error_G ; AEL_G] ;            CNR_error_G = [CNR_error_G ;ACNR_G] ;
            EL_error_B =[EL_error_B ; AEL_B] ;            CNR_error_B = [CNR_error_B ;ACNR_B] ;
        end     
    end
        %{
        %------------�Y�o��cycle slip�Ϸ|�����D�A�n�ⰻ��cyclslip���{������
        for r=1:AnGsat
            run_ADR(AGP_id(r),i) = Carrier_G( AGP_id(r) , 2 ) + ig(r , 1) - Tg(r , 1) + SCBg(r , 1) + relative_IGS(r , 1) - rec_bias ;
            run_trueRange(AGP_id(r),i) = ARr_G(r) ;
        end
        run_recbias(i)= rec_bias;
        for r=1:AnBsat
            run_ADR_B(ABD_id(r),i) = Carrier_B( ABD_id(r) , 2 ) + ib(r , 1) - Tb(r , 1) + SCBb(r , 1) - rec_bias_B + 0.225*i*lambda_B;
            run_trueRange_B(ABD_id(r),i) = ARr_B(r) ;
        end
        run_recbias_B(i)= rec_bias_B;
    end
    %}
    %----------------�ɶ�+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    %---------------------�p��T�����~�t(solid tide)----------------
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
    
    nsat_prior = nsat ; %   �����W�@��ìP����
    
    if first_P == 0 && all( Xxyz( i,: ) ) == 1  %�Ĥ@���w�즨�\��,������U�ɨ�
        first_P = i ;
    end
    nGsat_Mtrix( i ) = AnGsat ;
    nBsat_Mtrix( i ) = AnBsat ;
    nsat_Mtrix( i ) = Ansat ;
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

%---------------------��ƨ��ˮɶ�----------------------
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

Plotting_Rinex        % �p��STD, �e��

toc