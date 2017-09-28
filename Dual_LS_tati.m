tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List                   %   從Data_List.m內讀取資料
runtime = 7200 ;
Preparation               %   前置作業
Data  
obs_flag 

relay_station = runtime/2 ;
%---------------------初始設定--------------------------------
epk = 1 ;       %   累積data之秒數
ac_lim = 3       %   累積data之秒數的上限 ( 固定秒數為(ac_lim-1),浮動秒數為1 )
interval = 60        %   累積data間的間隔時刻數
first_P = 0 ;       %   用來記錄第一次成功定位的時刻數
first_PP = 0 ;
sat_var_SIGMA_all_ac = [] ;
CNR_Gac = [] ;           CNR_Bac = [] ;
Rho0_Gac = [] ;         Rho0_Bac = [] ;
EL_Gac = [] ;           EL_Bac = [] ;
AZ_Gac = [] ;           AZ_Bac = [] ;
sat_G_set = [] ;          sat_B_set = [] ;
SCBg_set = [] ;           SCBb_set = [] ;
relative_IGS_set = [] ;

bar = waitbar(0,'Please wait...');
for i = 1 : runtime
       
    str=['Positinoing ',num2str(i),'s'];    
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
                    SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
         
    wl_count = 0 ;                                      %計算GPS/BDS權重疊代次數
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
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    k = 1 ;                           %每一秒都把各顆BDS的資料依序抓進來
    
    while judge == 1
        
        while sat_EL( count ) < smaller_than_elevation
            count = count + 1 ;
        end
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
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
                
            elseif char( sat_sys( count ) ) == 66 %&& prn( count ) ~= 5 %&& ismember( prn(count) ,BD_sat ) == 1
               
                Rb( prn(count) , 1 ) = pr( count ) ;
                Rrb( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_B( prn(count) , 2 ) = ADR( count ) * lambda_B ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD衛星加上負號
                    BD_id_set( i,prn(count) ) = 1 ;
                    BD_id( k , 1) = prn( count ) ;
                    CNR_B( k , 1) = cnr0( count ) ;
                    Rho0_B( k , 1 ) = epr( count ) ;
                    Rho0_Bpr( k , 1 ) = pr( count ) ;
                    EL_B( k , 1 ) = sat_EL( count ) ;
                    AZ_B( k , 1 ) = sat_AZ( count ) ;
                    k = k + 1 ;
                elseif obs_flag == 3 && i > initial_time && count_B( prn(count),1 ) > converge_time && Carrier_B(prn(count),2) ~= 0 && Carrier_B(prn(count),1) ~= 0
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD衛星加上負號
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
      
    %----------------每次(每一秒)定位衛星的顆數---------------
    nGsat = length( GP_id ) ;    
    nBsat = length( BD_id ) ;
    nsat = nGsat + nBsat ;
    nGsat_Mtrix( i ) = nGsat ;
    nBsat_Mtrix( i ) = nBsat ;    
    nsat_Mtrix( i ) = nsat ; 
      
    %-----------------------------------------------------------------------------------------------    
    if nGsat ~= 0 && nBsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        gpst( i,: ) = r_gpst ;
                
        %----------------選擇觀測量------------------
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
        
        %-------------------納入current data----------------------------
        nGid( epk ) = nGsat ;
        nBid( epk ) = nBsat ;
        CNR_Gac = [ CNR_Gac ; CNR_G ] ;
        CNR_Bac = [ CNR_Bac ; CNR_B ] ;
        Rho0_Gac = [ Rho0_Gac ; Rho0_G ] ;
        Rho0_Bac = [ Rho0_Bac ; Rho0_B ] ;
        EL_Gac = [ EL_Gac ; EL_G ] ;
        EL_Bac = [ EL_Bac ; EL_B ] ;
        AZ_Gac = [ AZ_Gac ; AZ_G ] ;
        AZ_Bac = [ AZ_Bac ; AZ_B ] ;
        
        %-------------------計算covariance----------------------------
        EL = [ EL_Gac ; EL_Bac ] ;
        CNR = [ CNR_Gac ; CNR_Bac ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;        
        %-------------------選擇權重矩陣------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( ( sum(nGid)+sum(nBid) ) ) ;
        
        Rr_G = Rho0_G ;
        Rr_B = Rho0_B ;
        Delta_Mtrix( 1:3 ) = 1 ;
        
        while norm( Delta_Mtrix(1:3) ) > 0.01
            
            llh_0 = xyz2llh( rec_pos_0 ) ;
            lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
            lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
            height_rec = llh_0( 3 ) ;                  % unit : meter
    
            Tg = zeros( sum(nGid) , 1 ) ;
            ig = zeros( sum(nGid) , 1 ) ;
            Tb = zeros( sum(nBid) , 1 ) ;
            ib =  zeros( sum(nBid) , 1 ) ;
            
            %--------------------得到transmission time ( time_tx )----------------------
            [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
            [ XS_B , dtS_B , XS_tx_B , VS_tx_B , time_tx_B , no_eph_B , sys_idx_B ] = ...
                satellite_positions( r_gpst , Rr_B , BD_id , Eph_tatol_B , [] , [] , Tb , ib , rec_bias_satpos ) ;
            
            %------------利用 UNB3 模型計算對流層誤差-------------------------
            for n = 1 : sum(nGid)
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_Gac(n) * pi/180 ) ;
                Tg( n , 1 ) = RTROP ;                   
            end
            for n = 1 : sum(nBid)
                [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_Bac(n) * pi/180 ) ;
                Tb( n , 1 ) = RTROP_B ;
            end
            
            %------------使用GIM計算電離層誤差--------------------------
            for n = 1 : sum(nGid)
                [ iono , VTEC ]=...
                    GIM_corr( AZ_Gac(n) , EL_Gac(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                ig( n , 1 ) = iono ;
            end
            for n = 1 : sum(nBid)
                [ iono_B , VTEC_B ]=...
                    GIM_corr( AZ_Bac(n) , EL_Bac(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
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
            %----------------------使用saved data + current data---------------------------            
            sat_Gac = [ sat_G_set ; sat_G ] ;
            sat_Bac = [ sat_B_set ; sat_B ] ;
            
            %---------------------------修正地球自轉--------------------------------
            [ Xs_G , Rr_G_temp ] = Fix_Earth_Rotation( sat_Gac , rec_pos_0 ) ;
            [ Xs_B , Rr_B_temp ] = Fix_Earth_Rotation( sat_Bac , rec_pos_0 ) ;
            Rr_G = Rr_G_temp( end-nGsat+1:end ) ;
            Rr_B = Rr_B_temp( end-nBsat+1:end ) ;
            Xs = [ Xs_G ; Xs_B ] ;
            
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
            %----------------------使用saved data + current data---------------------------            
            SCBgac = [ SCBg_set ; SCBg ] ;
            SCBbac = [ SCBb_set ; SCBb ] ;
            relative_IGSac = [ relative_IGS_set ; relative_IGS ] ;
            
            %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------
            % 北斗的SCBb包含sat_bias,relative與group_delay
            Rhoc_G = Rho0_Gac - ig - Tg + SCBgac + relative_IGSac ;
            Rhoc_B = Rho0_Bac - ib - Tb + SCBbac ;
            Rhoc = [ Rhoc_G ; Rhoc_B ] ;
                       
            %-----------------------最小平方法定位-------------------------
            if nsat > 4
                
                [ Delta_Mtrix , Delta_Rho , V_hyber , G ] = LeastSquare_hyber_tati( Xs , rec_pos_0 , Rhoc , W , nGid , nBid ) ;
                
                while  ( abs(V_hyber(1)/V_hyber(2)) > G_B_condition || abs(V_hyber(2)/V_hyber(1)) > G_B_condition )                    
                    W = blkdiag( ( 1/V_hyber(1) )*W(1:sum(nGid) , 1:sum(nGid)) , ( 1/V_hyber(2) ) * W(sum(nGid)+1:end , sum(nGid)+1:end) ) ;
                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber , G ] = LeastSquare_hyber_tati( Xs , rec_pos_0 , Rhoc , W , nGid , nBid ) ;
                    wl_count = wl_count + 1 ;
                end
               
                rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                rec_bias_G = Delta_Mtrix(3+epk) ;
                rec_bias_B = Delta_Mtrix(end) ;
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
        rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置
                
        if first_P == 0 && all( Xxyz( i,: ) ) == 1  %第一次定位成功後,紀錄當下時刻
            first_P = i ;
        end
        cond = rem( i - first_P , interval ) ;
        
        %----------------------------倉庫矩陣(儲存舊時刻之資料)---------------------------------------
        if cond ~= 0            % 不符合區間時刻,清除current data
            CNR_Gac( end-nGsat+1:end ) = [] ;
            CNR_Bac( end-nBsat+1:end ) = [] ;
            Rho0_Gac( end-nGsat+1:end ) = [] ;
            Rho0_Bac( end-nBsat+1:end ) = [] ;
            EL_Gac( end-nGsat+1:end ) = [] ;
            EL_Bac( end-nBsat+1:end ) = [] ;
            AZ_Gac( end-nGsat+1:end ) = [] ;
            AZ_Bac( end-nBsat+1:end ) = [] ;            
        elseif  cond == 0 && epk < ac_lim            % 符合區間時刻且累積時刻量小於累積上限
            sat_G_set = sat_Gac ;
            sat_B_set = sat_Bac ;
            SCBg_set = SCBgac ;
            SCBb_set= SCBbac ;
            relative_IGS_set = relative_IGSac ;
            epk = epk + 1 ;
            first_P = i ;
        elseif  cond == 0 && epk == ac_lim            % 符合區間時刻且累積時刻量等於累積上限,清除oldest data
            CNR_Gac( 1:nGid(1) ) = [] ;
            CNR_Bac( 1:nBid(1) ) = [] ;
            Rho0_Gac( 1:nGid(1) ) = [] ;
            Rho0_Bac( 1:nBid(1) ) = [] ;
            EL_Gac( 1:nGid(1) ) = [] ;
            EL_Bac( 1:nBid(1) ) = [] ;
            AZ_Gac( 1:nGid(1) ) = [] ;
            AZ_Bac( 1:nBid(1) ) = [] ; 
            sat_G_set = [ sat_G_set( nGid(1)+1:end,: ) ; sat_G ] ;
            sat_B_set = [ sat_B_set( nBid(1)+1:end,: ) ; sat_B ] ;
            SCBg_set = [ SCBg_set( nGid(1)+1:end ) ; SCBg ] ;
            SCBb_set= [ SCBb_set( nBid(1)+1:end ) ; SCBb ] ;
            relative_IGS_set = [ relative_IGS_set( nGid(1)+1:end ) ; relative_IGS ] ;
            nGid(1) = [] ;
            nBid(1) = [] ;
            first_P = i ;
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
        
        if i > 1 && nsat ~= nsat_prior
            nsat_change(chg) = i ;
            chg = chg + 1 ;
        end
       
        nsat_prior = nsat ; %   紀錄上一秒衛星顆數
        
        if i == relay_station
            effTimeArray = sampling_time_start : relay_station ;
            rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;      %計算STD
            Xenu = zeros( 3600,3 ) ;
            for k = sampling_time_start : relay_station
                Xenu( k , : ) = xyz2enu( Xxyz( k , : ) , rec_pos_act )' ;
            end
            std = sqrt( var( Xenu(effTimeArray,:))  ) ;
            %-----------------四捨五入(round)取到小數點第X位---------------------
            std_enu = round( std .* 1e4 ) ./ 1e4 ;
            std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
            i_STD = [ relay_station ; std_enu' ; std_en' ]
            clear effTimeArray
        end
        
        if first_PP == 0 && all( Xxyz( i,: ) ) == 1  %第一次定位成功後,紀錄當下時刻
            first_PP = i ;
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
%sampling_time_start = 600 ;
sampling_time_end = runtime ;
effTimeArray = first_PP+1 : sampling_time_end ;

Plotting        % 計算STD, 畫圖

TaTi = [ ac_lim interval ]

toc