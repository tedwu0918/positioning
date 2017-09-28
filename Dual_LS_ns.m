tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List                   %   �qData_List.m��Ū�����
runtime = 7200 ;
Preparation               %   �e�m�@�~
Data  
obs_flag 
%---------------------��l�]�w--------------------------------
epk = 0 ;       %   �ֿndata�����,�q0�}�l
ac_lim = 30       %   �ֿndata����ƪ��W��
sat_var_SIGMA_all_ac = [] ;
CNR_Gac = [] ;           CNR_Bac = [] ;
Rho0_Gac = [] ;         Rho0_Bac = [] ;
EL_Gac = [] ;           EL_Bac = [] ;
AZ_Gac = [] ;           AZ_Bac = [] ;
sat_Gac = [] ;          sat_Bac = [] ;
SCBgac = [] ;           SCBbac = [] ;
relative_IGSac = [] ;

bar = waitbar(0,'Please wait...');
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];    
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %�C900��,���s��slagrange�����`�I
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
                    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX��(prn*data_num)���x�},�Y �Ĥ@�C���ìPprn1�b13�Ӯɨ誺x��m ,�ĤG�C���ìPprn2�b13�Ӯɨ誺x��m
         
    wl_count = 0 ;                                      %�p��GPS/BDS�v���|�N����
    id_temp_G = [] ;                                       id_temp_B = [] ;
    GP_id = [] ;                                                BD_id = [] ;
    Rho0_G = [] ;                                            Rho0_B = [] ;
    Rho0_Gpr = [] ;                                        Rho0_Bpr = [] ;
    CNR_G = [] ;                                             CNR_B = [] ;
    EL_G = [] ;                                                  EL_B = [] ;
    AZ_G = [] ;                                                 AZ_B = [] ;
    Rg = zeros( 32,1 ) ;                                    Rb = zeros( 14,1 ) ;
    Rrg( :,1 ) = Rrg( :,2 ) ;                                 Rrb( :,1 ) = Rrb( :,2 ) ;
    Rrg( :,2 ) = 0 ;                                             Rrb( :,2 ) = 0 ;
    Carrier_G( :,1 ) = Carrier_G( :,2 ) ;            Carrier_B( :,1 ) = Carrier_B( :,2 ) ;
    Carrier_G( :,2 ) = 0 ;                                  Carrier_B( :,2 ) = 0 ;
    Sat2s( :,1 ) = Sat2s( :,2 ) ;
    Sat2s( :,2 ) = 0 ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    k = 1 ;                           %�C�@����U��BDS����ƨ̧ǧ�i��
    
    while judge == 1
        
        while sat_EL( count ) < smaller_than_elevation
            count = count + 1 ;
        end
        
        if time( count ) == count_times                                                       % all( )==1 , �YAX�x�}�������������s
            if char( sat_sys( count ) ) == 80 &&  prn( count ) ~= 193 && all( AX( prn(count),: ) ) == 1 && AS( prn(count),1 ) ~= 999999.999999...
                    %&&  ismember( prn(count) ,GP_sat ) == 1
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_G( prn(count) , 2 ) = ADR( count ) * lambda_G ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    Sat2s( prn(count) , 2 ) = prn(count) ;
                    GP_id_set( i,prn(count) ) = 1 ;
                    GP_id( j , 1) = prn( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
                    Rho0_G( j , 1 ) = epr( count ) ;
                    Rho0_Gpr( j , 1 ) = pr( count ) ;
                    EL_G( j , 1 ) = sat_EL( count ) ;
                    AZ_G( j , 1 ) = sat_AZ( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                    Sat2s( prn(count) , 2 ) = prn(count) ;
                    GP_id_set( i,prn(count) ) = 1 ;
                    GP_id( j , 1) = prn( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
                    EL_G( j , 1 ) = sat_EL( count ) ;
                    AZ_G( j , 1 ) = sat_AZ( count ) ;
                    j = j+1 ;                   
                end
                
            elseif char( sat_sys( count ) ) == 66 %&&  prn( count ) ~= 5%&& ismember( prn(count) ,BD_sat ) == 1
                
                Rb( prn(count) , 1 ) = pr( count ) ;
                Rrb( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_B( prn(count) , 2 ) = ADR( count ) * lambda_B ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD�ìP�[�W�t��
                    BD_id_set( i,prn(count) ) = 1 ;
                    BD_id( k , 1) = prn( count ) ;
                    CNR_B( k , 1) = cnr0( count ) ;
                    Rho0_B( k , 1 ) = epr( count ) ;
                    Rho0_Bpr( k , 1 ) = pr( count ) ;
                    EL_B( k , 1 ) = sat_EL( count ) ;
                    AZ_B( k , 1 ) = sat_AZ( count ) ;
                    k = k + 1 ;
                elseif obs_flag == 3 && i > initial_time && count_B( prn(count),1 ) > converge_time && Carrier_B(prn(count),2) ~= 0 && Carrier_B(prn(count),1) ~= 0
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD�ìP�[�W�t��
                    BD_id_set( i,prn(count) ) = 1 ;
                    BD_id( k , 1) = prn( count ) ;
                    CNR_B( k , 1) = cnr0( count ) ;
                    EL_B( k , 1 ) = sat_EL( count ) ;
                    AZ_B( k , 1 ) = sat_AZ( count ) ;
                    k = k + 1 ;
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
    nGsat_Mtrix( i ) = nGsat ;
    nBsat_Mtrix( i ) = nBsat ;    
    nsat_Mtrix( i ) = nsat ; 
   
    %-----------------------------------------------------------------------------------------------    
    if nGsat ~= 0 && nBsat ~= 0
        
        %-----------�o�챵���ɶ�---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        gpst( i,: ) = r_gpst ;
                
        %----------------����[���q------------------
        if obs_flag == 1
            Rho0_G = Rho0_Gpr ;
            Rho0_B = Rho0_Bpr ;
        elseif obs_flag == 2
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
            [ Rho0_B , S_B , count_B ] = Doppler_Smoothed_Code( Rrb , Rb , S_B , BD_id , count_B , BMax ) ;
        elseif obs_flag == 3
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
            [ Rho0_B , S_B , count_B ] = Carrier_Smoothed_Code( Carrier_B , Rb , S_B , BD_id , count_B , BMax ) ;
        end
        
        if any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk < ac_lim
            epk = epk +1 ;      %   �Y�ìP���ƥ��ܰ�,�ֿn���+1
            CNR_Gac = [ CNR_Gac ; CNR_G ] ;
            CNR_Bac = [ CNR_Bac ; CNR_B ] ;
            Rho0_Gac = [ Rho0_Gac ; Rho0_G ] ;
            Rho0_Bac = [ Rho0_Bac ; Rho0_B ] ;
            EL_Gac = [ EL_Gac ; EL_G ] ;
            EL_Bac = [ EL_Bac ; EL_B ] ;
            AZ_Gac = [ AZ_Gac ; AZ_G ] ;
            AZ_Bac = [ AZ_Bac ; AZ_B ] ;
        elseif any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk == ac_lim
            CNR_Gac = [ CNR_Gac( nGsat+1:end ) ; CNR_G ] ;     %   �Y��F�ֿn���,�h�C�h�@���s�����,�N���ª��@����ƥ��
            CNR_Bac = [ CNR_Bac( nBsat+1:end ) ; CNR_B ] ;
            Rho0_Gac = [ Rho0_Gac( nGsat+1:end ) ; Rho0_G ] ;
            Rho0_Bac = [ Rho0_Bac( nBsat+1:end ) ; Rho0_B ] ;
            EL_Gac = [ EL_Gac( nGsat+1:end ) ; EL_G ] ;
            EL_Bac = [ EL_Bac( nBsat+1:end ) ; EL_B ] ;
            AZ_Gac = [ AZ_Gac( nGsat+1:end ) ; AZ_G ] ;
            AZ_Bac = [ AZ_Bac( nBsat+1:end ) ; AZ_B ] ;
        else
            epk = 1 ;       %   �Y�ìP�����ܰʫ�,���s�ֿn���
            CNR_Gac = CNR_G ;           CNR_Bac = CNR_B ;
            Rho0_Gac = Rho0_G ;         Rho0_Bac = Rho0_B ;
            EL_Gac = EL_G ;           EL_Bac = EL_B ;
            AZ_Gac = AZ_G ;           AZ_Bac = AZ_B ;
        end
               
        %-------------------�p��covariance----------------------------
        EL = [ EL_Gac ; EL_Bac ] ;
        CNR = [ CNR_Gac ; CNR_Bac ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        %-------------------����v���x�}------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( nsat*epk ) ;
        
        Rr_G = Rho0_G ;
        Rr_B = Rho0_B ;
        Delta_Mtrix( 1:3 ) = 1 ;
       
        while norm( Delta_Mtrix(1:3) ) > 0.01
            
            llh_0 = xyz2llh( rec_pos_0 ) ;
            lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
            lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
            height_rec = llh_0( 3 ) ;                  % unit : meter
    
            Tg = zeros( nGsat*epk , 1 ) ;
            ig = zeros( nGsat*epk , 1 ) ;
            Tb = zeros( nBsat*epk , 1 ) ;
            ib =  zeros( nBsat*epk , 1 ) ;
            
            %--------------------�o��transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
                satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
            [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
                satellite_positions( r_gpst , Rr_B , BD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
            
            %------------�Q�� UNB3 �ҫ��p���y�h�~�t-------------------------
            for n = 1 : nGsat*epk
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_Gac(n) * pi/180 ) ;
                Tg( n , 1 ) = RTROP ;                   
            end
            for n = 1 : nBsat*epk
                [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_Bac(n) * pi/180 ) ;
                Tb( n , 1 ) = RTROP_B ;
            end
            
            %------------�ϥ�GIM�p��q���h�~�t--------------------------
            for n = 1 : nGsat*epk
                [ iono , VTEC ]=...
                    GIM_corr( AZ_Gac(n) , EL_Gac(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                ig( n , 1 ) = iono ;
            end
            for n = 1 : nBsat*epk
                [ iono_B , VTEC_B ]=...
                    GIM_corr( AZ_Bac(n) , EL_Bac(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                ib( n , 1 ) = iono_B ;
            end
            
            %------------�ϥ� time_tx �p���K�P����m-----------------
            %--------------initialize----------------
            sat_G = zeros( nGsat , 3 ) ;
            SCBg = zeros( nGsat , 1 ) ;
            
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
                ephemeris_SCB = AS( GP_id(k),: ) ; %SCB=satellite clock bias
                SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
            end
            
            %-------------------�ϥ�stallie_position------------------------
            sat_B = XS_tx_B ;
            SCBb = c .* dtS_B ;
            
            if any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk ~= size(sat_G_set,1)/nGsat
                sat_Gac = [ sat_G_set ; sat_G ] ;          %   �p�G�ܮw�x�}����,�Y�ɨ�ƥ����W��,�h���ɨ褧��Ʒ|�ǤJ
                sat_Bac = [ sat_B_set ; sat_B ] ;
            elseif any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk == size(sat_G_set,1)/nGsat
                sat_Gac = [ sat_G_set( nGsat+1:end,: ) ; sat_G ] ;          %   �p�G�ܮw�x�}�w��,�Y�ɨ�Ƥw�F�W��,�h�|�����®ɨ誺���
                sat_Bac = [ sat_B_set( nBsat+1:end,: ) ; sat_B ] ;
            else
                sat_Gac = sat_G ;
                sat_Bac = sat_B ;
            end
            %---------------------------�ץ��a�y����--------------------------------
            [ Xs_G , Rr_G_temp ] = Fix_Earth_Rotation( sat_Gac , rec_pos_0 ) ;
            [ Xs_B , Rr_B_temp ] = Fix_Earth_Rotation( sat_Bac , rec_pos_0 ) ;
            Rr_G = Rr_G_temp( nGsat*(epk-1)+1:end ) ;
            Rr_B = Rr_B_temp( nBsat*(epk-1)+1:end ) ;
            Xs = [ Xs_G ; Xs_B ] ;
            
            %------------------�p��GPS�۹�׻~�t-------------------------------
            relative_IGS = zeros( nGsat , 1 ) ;
            tgd = 0 ;
            
            for n = 1 : nGsat
                icol = find_eph( Eph_tatol , GP_id( n ) , r_gpst ) ;
                Eph = Eph_tatol( : , icol ) ;
                
                [ satp , satv ] = satellite_orbits( time_tx( n ) , Eph , icol , [] ) ;
                
                Sp3_corr = -2*dot( satp , satv )/c ;
                relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
            end         
                                    
            if any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk ~= size(SCBg_set,1)/nGsat
                SCBgac = [ SCBg_set ; SCBg ] ;          %   �p�G�ܮw�x�}����,�Y�ɨ�ƥ����W��,�h���ɨ褧��Ʒ|�ǤJ
                SCBbac = [ SCBb_set ; SCBb ] ;
                relative_IGSac = [ relative_IGS_set ; relative_IGS ] ;
            elseif any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk == size(SCBg_set,1)/nGsat
                SCBgac = [ SCBg_set( nGsat+1:end ) ; SCBg ] ;          %   �p�G�ܮw�x�}�w��,�Y�ɨ�Ƥw�F�W��,�h�|�����®ɨ誺���
                SCBbac = [ SCBb_set( nBsat+1:end ) ; SCBb ] ;
                relative_IGSac = [ relative_IGS_set( nGsat+1:end ) ; relative_IGS ] ;
            else
                SCBgac = SCBg ;
                SCBbac = SCBb ;
                relative_IGSac = relative_IGS ;
            end       
            
            %--------------------�ץ������Z��(�U�����Ҭ�����)--------------------------
            % �_�檺SCBb�]�tsat_bias,relative�Pgroup_delay
            Rhoc_G = Rho0_Gac - ig - Tg + SCBgac + relative_IGSac ;
            Rhoc_B = Rho0_Bac - ib - Tb + SCBbac ;
            Rhoc = [ Rhoc_G ; Rhoc_B ] ;
                       
            %-----------------------�̤p����k�w��-------------------------
            if nsat > 4
                
                [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber_ns( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , epk ) ;
                
                while  ( abs(V_hyber(1)/V_hyber(2)) > G_B_condition || abs(V_hyber(2)/V_hyber(1)) > G_B_condition )                    
                    W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat*epk , 1:nGsat*epk) , ( 1/V_hyber(2) ) * W(nGsat*epk+1:end , nGsat*epk+1:end) ) ;                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber_ns( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , epk ) ;
                    wl_count = wl_count + 1 ;    
                end
                
                rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                rec_bias_G = rec_bias_G + Delta_Mtrix(3+epk) ;
                rec_bias_B = rec_bias_B + Delta_Mtrix(3+2*epk) ;
                rec_pos_0 = rec_pos ;
            else
                break
            end       
            
        end
        
        if nsat < 5
            Xxyz( i , : ) = rec_pos_0 ;
        else
            Xxyz( i , : ) = rec_pos ;
        end        
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   �N�ұo����m��s���U�@�ɨ褧��l��m
            
            %----------------------------�ܮw�x�}(�x�s�®ɨ褧���)---------------------------------------
        if any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk ~= size(sat_G_set,1)/nGsat
            sat_G_set = [ sat_G_set ; sat_G ] ;          %   �p�G�ܮw�x�}����,�Y�ɨ�ƥ����W��,�h���ɨ褧��Ʒ|�ǤJ
            sat_B_set = [ sat_B_set ; sat_B ] ;
            SCBg_set = [ SCBg_set ; SCBg ] ;
            SCBb_set= [ SCBb_set ; SCBb ] ;
            relative_IGS_set = [ relative_IGS_set ; relative_IGS ] ;
        elseif any( Sat2s(:,2) - Sat2s(:,1) ) == 0 && epk == size(sat_G_set,1)/nGsat
            sat_G_set = [ sat_G_set( nGsat+1:end,: ) ; sat_G ] ;          %   �p�G�ܮw�x�}�w��,�Y�ɨ�Ƥw�F�W��,�h�|�����®ɨ誺���
            sat_B_set = [ sat_B_set( nBsat+1:end,: ) ; sat_B ] ;
            SCBg_set = [ SCBg_set( nGsat+1:end ) ; SCBg ] ;
            SCBb_set= [ SCBb_set( nBsat+1:end ) ; SCBb ] ;
            relative_IGS_set = [ relative_IGS_set( nGsat+1:end ) ; relative_IGS ] ;
        else
            sat_G_set = sat_G ;
            sat_B_set = sat_B ;
            SCBg_set = SCBg ;
            SCBb_set = SCBb ;
            relative_IGS_set = relative_IGS ;
        end
            
    end    
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
%sampling_time_start = 600 ;
sampling_time_end = runtime ;
effTimeArray = sampling_time_start : sampling_time_end ;

Plotting        % �p��STD, �e��

toc