tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List           %   從Data_List.m內讀取資料
runtime = 10000 ;
Preparation               %   前置作業

T = 1; % positioning interval
% State vector is as [x Vx y Vy z Vz b d].', i.e. the coordinate (x,y,z),
% the clock bias b, and their derivatives.

% Kalman Parameters initial setting
sigma= 1 ;      Rhoerror = 4 ;       Pvalue = 10 ;
%Sb = (1.1e-19)*c^2 ;        Sd = (4.3e-20)*c^2 ;     
%Sb = (4e-19)*c^2 ;        Sd = ((pi^2)*16e-20)*c^2 ;%The study of GPS Time Transfer Based on Extended Kalman Filter
Sb = (2e-19)*c^2 ;          Sd = ((pi^2)*56e-20)*c^2 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)
%Sb = (4e-19)*c^2 ;        Sd = (1.58e-18)*c^2 ;    

Xp = zeros( 10 , 1 ) ;
Pp = eye( 10 )*Pvalue ;
Xu = zeros( 10 , 1 ) ;

Qc = [Sb*T+Sd*T*T*T/3     Sd*T*T/2 ;
                    Sd*T*T/2                Sd*T ] ;
Qxyz = sigma^2 * [T^3/3      T^2/2 ;
                                    T^2/2        T      ] ;
Q = blkdiag( Qxyz , Qxyz , Qxyz , Qc , Qc ) ;

ff = [ 1 T ; 0 1 ] ;
fy = blkdiag( ff , ff , ff , ff , ff ) ;

first_P = 1 ;
KF = 0 ;
%-------------------------------------------------------------------------------------
bar = waitbar(0,'Please wait...');   
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];    
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
    
    GP_id = [] ;                        BD_id = [] ;
    Rho0_G = [] ;                    Rho0_B = [] ;
    Rho0_Gpr = [] ;                Rho0_Bpr = [] ;
    CNR_G = [] ;                     CNR_B = [] ;
    EL_G = [] ;                         EL_B = [] ;
    AZ_G = [] ;                        AZ_B = [] ;
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
                    %&& prn(count)~=31 && prn(count)~=21 && prn(count)~=20
                
            Rg( prn(count) , 1 ) = pr( count ) ;
            Rrg( prn(count) , 2 ) = pr_rate( count ) ;
            Carrier_G( prn(count) , 2 ) = ADR( count ) * lambda_G ;
            ADR_G( i,prn(count) ) = ADR( count ) ;
            
            if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                GP_id_set( i,prn(count) ) = 1 ;
                GP_id( j , 1) = prn( count ) ;
                CNR_G( j , 1) = cnr0( count ) ;
                Rho0_G( j , 1 ) = epr( count ) ;
                Rho0_Gpr( j , 1 ) = pr( count ) ;
                EL_G( j , 1 ) = sat_EL( count ) ;
                AZ_G( j , 1 ) = sat_AZ( count ) ;
                j = j+1 ;
            elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                GP_id_set( i,prn(count) ) = 1 ;
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
                ADR_B( i,prn(count) ) = ADR( count ) ;
                
                if obs_flag <= 2 || ( obs_flag == 3 && i <= initial_time )
                    BD_id_set( i,prn(count) ) = 1 ;
                    BD_id( k , 1) = prn( count ) ;
                    CNR_B( k , 1) = cnr0( count ) ;
                    Rho0_B( k , 1 ) = epr( count ) ;
                    Rho0_Bpr( k , 1 ) = pr( count ) ;
                    EL_B( k , 1 ) = sat_EL( count ) ;
                    AZ_B( k , 1 ) = sat_AZ( count ) ;
                    k = k + 1 ;
                elseif obs_flag == 3 && i > initial_time && count_B( prn(count),1 ) > converge_time && Carrier_B(prn(count),2) ~= 0 && Carrier_B(prn(count),1) ~= 0
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
    
    %--------------------------------------------------------------------    
    if nGsat ~= 0 && nBsat ~= 0
        
       %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %-------------------計算covariance----------------------------
        EL = [ EL_G ; EL_B ] ;
        CNR = [ CNR_G ; CNR_B ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;   
        
        %-------------------選擇權重矩陣------------------------------   
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        
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
        
        Rr_G = Rho0_G ;
        Rr_B = Rho0_B ;
        Delta_Mtrix( 1:5 ) = 1 ;
        
        while norm( Delta_Mtrix(1:3) ) > 1e-2
            
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
                [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_G(n) * pi/180 ) ;
                Tg( n , 1 ) = RTROP ;                   % lat_rec : radians
            end
            for n = 1 : nBsat
                [ RTROP_B , HZD_B , HMF_B , WZD_B , WMF_B ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , EL_B(n) * pi/180 ) ;
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
                sat_x_temp = LagrangeInter( Data_time , AX( GP_id(k),: ) , time_tx( k ) ) ;
                sat_y_temp = LagrangeInter( Data_time , AY( GP_id(k),: ) , time_tx( k ) ) ;
                sat_z_temp = LagrangeInter( Data_time , AZ( GP_id(k),: ) , time_tx( k ) ) ;
                sat_G( k , 1 ) = sat_x_temp ;
                sat_G( k , 2 ) = sat_y_temp ;
                sat_G( k , 3 ) = sat_z_temp ;
                for t = 1 : length( AX( GP_id(k),: ) )
                    if Data_time( t ) - time_tx( k ) > 0
                        Interval = t-1 : t ;
                        break
                    end
                end
                ephemeris_SCB = AS( GP_id(k),: ) ; %SCB=satellite clock bias
                SCB_temp = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
                SCBg( k , 1 ) = SCB_temp ;
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
            if EPR_flag == 1
                Rhoc_G = Rho0_G + SCBg + relative_IGS ;
                Rhoc_B = Rho0_B + SCBb ;
            else
                Rhoc_G = Rho0_G - ig - Tg + SCBg + relative_IGS ;
                Rhoc_B = Rho0_B - ib - Tb + SCBb ;
            end
            Rhoc = [ Rhoc_G ; Rhoc_B ] ;
            
            %---------------------------修正地球自轉--------------------------------
            [ Xs_G , Rr_G ] = Fix_Earth_Rotation( sat_G , rec_pos_0 ) ;
            [ Xs_B , Rr_B ] = Fix_Earth_Rotation( sat_B , rec_pos_0 ) ;
            Xs = [ Xs_G ; Xs_B ] ;
            
            if nsat > 4
                
                [ Delta_Mtrix , Delta_Rho , V_hyber  ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                
                while  ( abs( V_hyber(1)/V_hyber(2) )  > G_B_condition ||  abs( V_hyber(2)/V_hyber(1) ) > G_B_condition )
                    W = blkdiag( ( 1/V_hyber(1) )*W(1:nGsat , 1:nGsat) , ( 1/V_hyber(2) ) * W(nGsat+1:nsat , nGsat+1:nsat) ) ;
                    
                    [ Delta_Mtrix , Delta_Rho , V_hyber ] = LeastSquare_hyber( Xs , rec_pos_0 , Rhoc , W , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
                    wl_count = wl_count + 1 ;
                end
                
                rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                rec_bias_G = rec_bias_G + Delta_Mtrix(4) ;
                rec_bias_B = rec_bias_B + Delta_Mtrix(5) ;
                rec_pos_0 = rec_pos ;
            else
                break
            end            
        end
        
        if first_P == 1 && nsat < 5
            Xxyz( i , : ) = rec_pos_0 ;
        elseif first_P == 1
            Xxyz( i , : ) = rec_pos ;  
            Xp([1,3,5]) = Delta_Mtrix(1:3).' ;
            Xp(7) = Delta_Mtrix(4) ;
            Xp(9) = Delta_Mtrix(5) ;
            first_P = 0 ;
            KF = 1 ;
        end
        
        %--------------------Kalman Filter------------------------------
        if KF == 1 %&& i > 1 && any( Xxyz(i-1,:) ) == 1
            R = eye( size( Xs , 1 ) ) * Rhoerror;           
            [ Val , H ] = H_Mtrix( rec_pos_0 , Xs , nGsat , nBsat , rec_bias_G , rec_bias_B ) ;
            Z = Rhoc - Val ;    % Z= z - zp
            
            K = Pp * H' / (H * Pp * H.' + R) ;
            
            Xu = Xp + K * ( Z - H*Xp ) ;        %   Z為修正過後的pseudorange
            I = eye( size(Xu,1),size(Xu,1) ) ;
            Pu = (I - K * H) * Pp;
            
            Xp = fy * Xu ;                  % fy為State transition matrix
            Pp = fy * Pu * fy.' + Q ;
            
            rec_pos = rec_pos_0 + Xu( [1,3,5] )' ;
            rec_bias_G = rec_bias_G + Xu(7) ;
            rec_bias_B = rec_bias_B + Xu(9) ;
            Xxyz( i , : ) = rec_pos ;
            
        end        
        
        %---------------count DOP--------------------
        [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
        rec_pos_0 = Xxyz( i , : ) ;             %   將所得之位置更新為下一時刻之初始位置
    end
    
        %----------------時間+1-------------------------
        time_week_last_epoch = time_week ;
        time_week = time( count ) ;
        
        time_day_UTC = time_day_UTC + 1 ;
        local_time = local_time + 1 ;
        
        %---------------------計算固體潮汐誤差(solid tide)----------------
        if Data ~= no_solid_tide
            solid_tide_error( i , 1 ) = interp1( time_solid , e_solid , time_day_UTC ) ;
            solid_tide_error( i , 2 ) = interp1( time_solid , n_solid , time_day_UTC ) ;
            solid_tide_error( i , 3 ) = interp1( time_solid , u_solid , time_day_UTC ) ;
        end
        
        if i > 1 && nsat ~= nsat_prior
            nsat_change(chg) = i ;
            chg = chg + 1 ;
        end
        
        nsat_prior = nsat ;     %   紀錄上一秒衛星顆數
        
end
close( bar ) ;

switch obs_flag
    case { 0 }
        obs = 'Dual_EPR_errKF' ;
    case { 1 }
        obs = 'Dual_PR_errKF' ;
    case { 2 }
        obs = 'Dual_DSC_errKF' ;
    case { 3 }
        obs = 'Dual_CSC_errKF' ;
end

%---------------------資料取樣時間----------------------
%sampling_time_start = 100 ;
sampling_time_end = runtime ;
effTimeArray = sampling_time_start : sampling_time_end ;

%---------------------顯示參數設定----------------------
[ Data_start_end_obs_elevation_Q_R ] = { Data sampling_time_start sampling_time_end obs ...
    smaller_than_elevation sigma Rhoerror }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;      %計算STD

Xenu = zeros( runtime , 3 ) ;
for i = 1 : runtime 
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( var( Xenu(effTimeArray,:))  ) ;
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ; 
STD = [ std_enu' ; std_en' ] 

%rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
%rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ; 
%RMS = [ rms_enu' ; rms_en' ] ;

%error = [ std' , rms' ] ;
%Error = [ STD , RMS ]

%--------------------2D plot----------------------
figure( 10 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
title( 'Scatter Plot (GPS/BDS)' ) ;
print( '-dpng',  'Scatter plot_KF' , '-r600'  ) ;       %Change "-r600" to the required DPI
%{
%-------------------3D plot------------------------------
    figure( 11 );
    plot(effTimeArray , Xenu( effTimeArray , 1 ) , 'r-.' ,...
            effTimeArray , Xenu( effTimeArray , 2 ) , 'b-.' ,...
            effTimeArray , Xenu( effTimeArray , 3 ) , 'g-.' ) ;
    grid on;
    xlabel('time(s)');              
    ylabel('Errors(m)');
    legend('East','North','UP');
    print( '-dpng',  'ENU error' , '-r600'  ) ; 
%}
    figure( 12 );
    subplot(3,1,1);
    plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
    xlabel('Time (s)');                 
    ylabel('EAST (m)');
    title('ENU 位置與時間關係');
    grid on;                        hold on;

    subplot(3,1,2);
    plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
    xlabel('Time (s)');                 
    ylabel('NORTH (m)');
    grid on;                        hold on;

    subplot(3,1,3);
    plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
    xlabel('Time (s)');                 
    ylabel('UP (m)');
    grid on;                        hold on;
    print( '-dpng',  'ENU_KF' , '-r600' ) ;       %Change "-r600" to the required DPI
    toc