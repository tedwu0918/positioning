tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List                   %   從Data_List.m內讀取資料
runtime = 7200 ;   % ntou0322 : 47025
Preparation              %   前置作業
Data
for_k = 1 
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
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_G( prn(count) , 2 ) = ADR( count ) * lambda_G ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    Sat2s( prn(count) , 2 ) = prn(count) ;
                    GP_id( j , 1) = prn( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
                    Rho0_G( j , 1 ) = epr( count ) ;
                    Rho0_Gpr( j , 1 ) = pr( count ) ;
                    EL_G( j , 1 ) = sat_EL( count ) ;
                    AZ_G( j , 1 ) = sat_AZ( count ) ;
                    j = j+1 ;
                elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                    Sat2s( prn(count) , 2 ) = prn(count) ;
                    GP_id( j , 1) = prn( count ) ;
                    CNR_G( j , 1) = cnr0( count ) ;
                    EL_G( j , 1 ) = sat_EL( count ) ;
                    AZ_G( j , 1 ) = sat_AZ( count ) ;
                    j = j+1 ;
                end
                
            elseif char( sat_sys( count ) ) == 66 && prn( count ) ~= 5
                
                Rb( prn(count) , 1 ) = pr( count ) ;
                Rrb( prn(count) , 2 ) = pr_rate( count ) ;
                Carrier_B( prn(count) , 2 ) = ADR( count ) * lambda_B ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD衛星加100
                    BD_id( k , 1) = prn( count ) ;
                    CNR_B( k , 1) = cnr0( count ) ;
                    Rho0_B( k , 1 ) = epr( count ) ;
                    Rho0_Bpr( k , 1 ) = pr( count ) ;
                    EL_B( k , 1 ) = sat_EL( count ) ;
                    AZ_B( k , 1 ) = sat_AZ( count ) ;
                    k = k + 1 ;
                elseif obs_flag == 3 && i > initial_time && count_B( prn(count),1 ) > converge_time && Carrier_B(prn(count),2) ~= 0 && Carrier_B(prn(count),1) ~= 0
                    Sat2s( prn(count)+32 , 2 ) = prn(count) + 100 ;      %   BD衛星加100
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
    
    if i > 1 && nsat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    nsat_prior = nsat ;     %   紀錄上一秒衛星顆數
    %-----------------------------------------------------------------------------------------------    
    if nGsat ~= 0 && nBsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        gpst( i,: ) = r_gpst ;
        gpsTime = r_gpst ;
        
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
                
        Temp = [ cat(1,GP_id,BD_id+100) , cat(1,CNR_G,CNR_B) , cat(1,Rho0_G,Rho0_B) , cat(1,EL_G,EL_B) , cat(1,AZ_G,AZ_B) ] ;  %   當下時刻之資料
  
        %-------------------------------------由衛星顆數判斷PBF之case-----------------------------------------
        if i > 2 && any( Xxyz(i-1,:) ) == 1 && any( Xxyz(i-2,:) ) == 1 && any( Sat2s(:,2) - Sat2s(:,1) ) ~= 0 %case:衛星顆數變動(any:判斷矩陣內任一元素是否為非零,只要有非零值,傳回1)
            clear gpsTime Taf Tbf TcPri Tc D gain loss want
            if all( ( Sat2s(:,2) - Sat2s(:,1) )>=0 ) == 1         %   case:衛星顆數變多 ( all : 判斷矩陣內是否全為非零值,全部非零,傳回1 )                
                recPos = [ Xxyz(i-1,:) ; Xxyz(i-1,:) ] ;               % 1. before(t)   2. after(t)
                recBias = [ RecBias(i-1,1) RecBias(i-1,2) ; RecBias(i-1,1) RecBias(i-1,2) ] ;
                gpsTime = [ gpst(i) , gpst(i) ] ;
                gain = Sat2s( find(( Sat2s(:,2) - Sat2s(:,1) ) ~=0) , 2 ) ;
                [ output , index ] = ismember( gain , Temp(:,1) ) ;
                want = setdiff( 1:size(Temp(:,1)), index ) ;            %   c = setdiff(a,b) , 返回屬於a但不屬於b的不同元素的集合，C = a-b
                Tbf = Temp( want,: ) ;
                D = { Tbf(:,1) , Tbf(:,2) ,Tbf(:,3) , Tbf(:,4) , Tbf(:,5) ; Temp(:,1) , Temp(:,2) , Temp(:,3) , Temp(:,4) , Temp(:,5) } ;%當下時刻之資料最後跑
                Gain( 1,G_count ) = i ;     %   紀錄"衛星顆數增加"出現幾次
                G_count = G_count + 1 ;
            elseif all( ( Sat2s(:,2) - Sat2s(:,1) )<=0 ) == 1         %   case:衛星顆數變少 ( all : 判斷矩陣內是否全為非零值,全部非零,傳回1 )
                recPos = [ Xxyz(i-2,:) ; Xxyz(i-2,:) ; Xxyz(i-1,:) ] ;  % 1. before(t-1)   2. after(t-1)    3. after(t)
                recBias = [ RecBias(i-2,1) RecBias(i-2,2) ; RecBias(i-2,1) RecBias(i-2,2) ; RecBias(i-1,1) RecBias(i-1,2) ] ;
                gpsTime = [ gpst(i-1) , gpst(i-1) , gpst(i) ] ;
                loss = Sat2s( find(( Sat2s(:,2) - Sat2s(:,1) ) ~=0) , 1 ) ;
                [ output , index ] = ismember( loss , TempPri(:,1) ) ;
                want = setdiff( 1:size(TempPri(:,1)), index ) ;
                Taf = TempPri( want,: ) ;
                D = { Taf(:,1) , Taf(:,2) ,Taf(:,3) , Taf(:,4) , Taf(:,5) ; TempPri(:,1) , TempPri(:,2) , TempPri(:,3) , TempPri(:,4) , TempPri(:,5) ;...
                    Temp(:,1) , Temp(:,2) , Temp(:,3) , Temp(:,4) , Temp(:,5) } ;%上一秒after, 上一秒之before, 當秒資料
                Loss( 1,L_count ) = i  ;     %   紀錄"衛星顆數減少"出現幾次
                L_count = L_count + 1 ;
            else    %    case:衛星顆數有增有減                
                recPos = [ Xxyz(i-2,:) ; Xxyz(i-2,:) ; Xxyz(i-1,:) ; Xxyz(i-1,:) ] ;% 1. commo(t-1)   2. before(t-1)  3. common(t)   4.after(t)-->( 當下資料 )
                recBias = [ RecBias(i-2,1) RecBias(i-2,2) ; RecBias(i-2,1) RecBias(i-2,2) ; RecBias(i-1,1) RecBias(i-1,2) ; RecBias(i-1,1) RecBias(i-1,2) ] ;
                gpsTime = [ gpst(i-1) , gpst(i-1) , gpst(i) , gpst(i) ] ;
                com = intersect( Sat2s(:,2) , Sat2s(:,1) ) ;     %   找出common的集合
                common = com( find( com ~= 0 ) ) ;
                [ output , index ] = ismember( common , TempPri(:,1) ) ;      %   找出 t-1 秒時的common的index
                want = sort( index ) ;
                TcPri = TempPri( want,: ) ;    %   將 t-1 秒common的資料存成TcPri
                [ output , index ] = ismember( common , Temp(:,1) ) ;      %   找出 t 秒時的common的index
                want = sort( index ) ;
                Tc = Temp( want,: ) ;    %   將 t 秒common的資料存成Tc
                Tbf = TempPri ;    %   將 t-1 秒before的資料存成Tbf
                D = { TcPri(:,1) , TcPri(:,2) ,TcPri(:,3) , TcPri(:,4) , TcPri(:,5) ; Tbf(:,1) , Tbf(:,2) , Tbf(:,3) , Tbf(:,4) , Tbf(:,5) ; Tc(:,1) , Tc(:,2) , Tc(:,3) , Tc(:,4) , Tc(:,5) ; ...
                    Temp(:,1) , Temp(:,2) , Temp(:,3) , Temp(:,4) , Temp(:,5) } ;%上一秒common, 上一秒之before, 當秒common, 當秒after(當下資料)
                GainLoss( 1,GL_count ) = i ;     %   紀錄"衛星顆數有增有減"出現幾次
                GL_count = GL_count + 1 ;
            end
        end
        
        for ii = 1 : size(gpsTime,2)
            clear Delta
            del = 1 ;
            if size(gpsTime,2) > 1
                clear iG GP_id CNR_G Rho0_G EL_G AZ_G
                clear iB BD_id CNR_B Rho0_B EL_B AZ_B
                iG = find( D{ ii,1 }<100 ) ;
                iB = find( D{ ii,1 }>100 ) ;
                GP_id = D{ ii,1 }(iG);
                BD_id = D{ ii,1 }(iB) - 100 ;
                CNR_G = D{ ii,2 }(iG) ;
                CNR_B = D{ ii,2 }(iB) ;
                Rho0_G = D{ ii,3 }(iG) ;
                Rho0_B = D{ ii,3 }(iB) ;
                EL_G = D{ ii,4 }(iG) ;
                EL_B = D{ ii,4 }(iB) ;
                AZ_G = D{ ii,5 }(iG) ;
                AZ_B = D{ ii,5 }(iB) ;
                rec_pos_0 = recPos( ii,: ) ;
                r_gpst = gpsTime( ii ) ;
                rec_bias_G = recBias( ii,1 ) ;
                rec_bias_B = recBias( ii,2 ) ;
                
                %----------------每次(每一秒)定位衛星的顆數---------------
                nGsat = length( GP_id ) ;
                nBsat = length( BD_id ) ;
                nsat = nGsat + nBsat ;
            end
            
            %-------------------計算covariance----------------------------
            EL = [ EL_G ; EL_B ] ;
            CNR = [ CNR_G ; CNR_B ] ;
            [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
            
            %-------------------選擇權重矩陣------------------------------
            W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
            %W = eye( nsat ) ;
            
            Rr_G = Rho0_G ;
            Rr_B = Rho0_B ;
            Delta_Mtrix( 1:3 ) = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 0.01
                
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
                    
                    [ Delta_Mtrix , Delta_Rho , G , V_hyber  ] = LeastSquare_hyber_PBF( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                    
                    while  ( abs(V_hyber(1)/V_hyber(2)) > G_B_condition || abs(V_hyber(2)/V_hyber(1)) > G_B_condition ) 
                        W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat , 1:nGsat) , ( 1/V_hyber(2) ) * W(nGsat+1:nsat , nGsat+1:nsat) ) ;
                        
                        [ Delta_Mtrix , Delta_Rho , G , V_hyber  ] = LeastSquare_hyber_PBF( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                        wl_count = wl_count + 1 ;
                    end
                    rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias_G = rec_bias_G + Delta_Mtrix(4) ;
                    rec_bias_B = rec_bias_B + Delta_Mtrix(5) ;
                    rec_pos_0 = rec_pos ;
                    
                    Delta( :,del ) = Delta_Mtrix ;
                    del = del + 1 ;
                    %Delta_set(i,:) = Delta_set(i,:) + Delta_Mtrix' ;
                    
                else
                    break
                end
                
            end
            Delta_total( i,ii ) = { Delta } ;
        end
        
        if nsat < 5
            Xxyz( i , : ) = rec_pos_0 ;
        else
            Xxyz( i , : ) = rec_pos ;
        end    
        
        %----------------------------Position Bias Filter ( PBF ) ----------------------------------
        dx = zeros( 3,1 ) ;
        if ii == 2
            dx = sum( Delta_total{i,2}(1:3,:) , 2 ) - sum( Delta_total{i,1}(1:3,:) , 2 ) ;  %   after - before   
        elseif ii == 3      
            dx = sum( Delta_total{i,1}(1:3,:) , 2 ) - sum( Delta_total{i,2}(1:3,:) , 2 ) ;  %   after - before
        elseif ii == 4
            dxGain = sum( Delta_total{i,4}(1:3,:) , 2 ) - sum( Delta_total{i,3}(1:3,:) , 2 ) ;    %   gain : after(t) - common(t)
            dxLoss = sum( Delta_total{i,1}(1:3,:) , 2 ) - sum( Delta_total{i,2}(1:3,:) , 2 ) ;    %   loss : common(t-1) - before(t-1)
            dx = dxGain + dxLoss ;            
        end
        Cx = for_k * Cx + dx' ;
        Xxyz( i , : ) = Xxyz( i , : ) - Cx ;
        Cx_set( i,: ) = Cx ;
        
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;        
        rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置
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
    
    RecBias( i,1 ) = rec_bias_G ;   %   存下每一秒的接收機時錶誤差
    RecBias( i,2 ) = rec_bias_B ;
    
    if nGsat ~= 0 && nBsat ~= 0
        TempPri = Temp ;
        clear Temp
    end
    
end

close( bar ) ;

switch obs_flag
    case { 0 }
        obs = 'EPR-PBF' ;
    case { 1 }
        obs = 'PR-PBF' ;
    case { 2 }
        obs = 'DSC-PBF' ;
    case { 3 }
        obs = 'CSC-PBF' ;
end
%---------------------資料取樣時間----------------------
sampling_time_end = runtime ;
effTimeArray = sampling_time_start : sampling_time_end ;

for_k
Plotting        % 計算STD, 畫圖

ttt = 1 : i  ;

figure( 13 );
subplot(3,1,1);
plot( ttt , Cx_set( ttt , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('X (m)');
title('Accumulated Position Bias Vector ( Cx )');
grid on;                        hold on;

subplot(3,1,2);
plot(ttt, Cx_set( ttt , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('Y (m)');
grid on;                        hold on;

subplot(3,1,3);
plot(ttt, Cx_set( ttt , 3 ) , 'g-' ) ;
xlabel('Time (s)');
ylabel('Z (m)');
grid on;                        hold on;
print( '-dpng',  'Cx_plot' , '-r600' ) ;       %Change "-r600" to the required DPI

toc