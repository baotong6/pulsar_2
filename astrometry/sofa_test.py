# -*- coding: utf-8 -*-
import numpy as np

path='/Users/baotong/desktop/pulsar/'
correction3=np.genfromtxt(path+'correction3.txt')
def r_CIS(r_CTS,JD_TT):
    def de_ra(x):
        return x*np.pi/180.0
    #IAU2000B章动系数
    # -------------------------------------------------------------------------------------------
    #常值参数
    #delta_AT = 35; #协调世界时换算原子时的改正值
    delta_T = 68.2; #世界时换算力学时的改正值
    xp = 0.14086;   #瞬时极在地级坐标系中的坐标/度
    yp = -0.676091; #瞬时极在地级坐标系中的坐标/度
    omega_e = [0,0,7.292115e-5];  #地球自转角速度
    # -------------------------------------------------------------------------------------------
    #时间计算
    #day_TAI = day_UTC + delta_AT/86400;#国际原子时TAI
    #地球时TT
    #day_TT  = day_TAI + 32.184/86400;   #地球时TT
    #JD_TT =2457061.99
    T_TT = (JD_TT-2451545.0)/36525;     #自J2000TT起算的儒略世纪数

    #质心动力学时TDB
    JD_TDB = JD_TT + 0.001657*np.sin(de_ra(628.3076*T_TT+6.2401))+0.000022*np.sin(de_ra(575.3385*T_TT+4.2970))\
              +0.000014*np.sin(de_ra(1256.6152*T_TT+6.1969)) \
              +0.000005*np.sin(de_ra(606.9777*T_TT+4.0212))\
              +0.000005*np.sin(de_ra(52.9691*T_TT+0.4444))\
              +0.000002*np.sin(de_ra(21.3299*T_TT+5.5431))\
              +0.000010*T_TT*np.sin(de_ra(628.3076*T_TT+4.2490));

    T_TDB = (JD_TDB-2451545.0)/36525;   #自J2000TDB起算的儒略世纪数
    #世界时UT1
    JD_UT1 = JD_TT - delta_T/86400;   #世界时UT1
    DU = JD_UT1 - 2451545.0;            #自2000年1月1日12h UT1起算的UT1天数

    #---------------------------------------------------------------------------------------------
    def de_ra(x):
        return x*np.pi/180.0
    def de_ra(x):
        return x*np.pi/180.0

    def M1(x):
        x=de_ra(x)
        fun=np.array([[1,0,0],[0,np.cos(x),np.sin(x)],[0,-np.sin(x),np.cos(x)]])
        return fun

    def M2(x):
        x = de_ra(x)
        fun = np.array([[np.cos(x),0,-np.sin(x)],[0,1,0],[np.sin(x),0,np.cos(x)]])
        return fun

    def M3(x):
        x = de_ra(x)
        fun = np.array([[np.cos(x),np.sin(x),0],[-np.sin(x),np.cos(x),0],[0,0,1]])
        return fun

    def CIS_MOD(T_TDB):
        zetaA  = 2.650545 + 2306.083227*T_TDB + 0.2988499*T_TDB**2 + 0.01801828*T_TDB**3- 0.000005971*T_TDB**4 - 0.0000003173*T_TDB**5;
        thetaA = 2004.191903*T_TDB - 0.4294934*T_TDB**2 - 0.04182264*T_TDB**3- 0.000007089*T_TDB**4 - 0.0000001274*T_TDB**5;
        zA = -2.650545 + 2306.077181*T_TDB + 1.0927348*T_TDB**2 + 0.01826837*T_TDB**3- 0.000028596*T_TDB**4 - 0.0000002904*T_TDB**5;

        CIS2MOD = M3(-zA/3600).dot((M2(thetaA/3600)).dot(M3(-zetaA/3600)));
        #协议天球坐标系to瞬时平天球坐标系
        return CIS2MOD

    def MOD_TOD(T_TDB):
        delta_posai   = -0.135e-3*1e7;
        delta_epsilon = 0.388e-3*1e7;
        phi = np.array([[485868.249036,1717915923.2178,31.8792,0.051635,-0.00024470],
           [1287104.79305,129596581.0481,-0.5532,0.000136,-0.00001149],
           [335779.526232,1739527262.8478,-12.7512,-0.001037,0.00000417],
           [1072260.70369,1602961601.2090,-6.3706,0.006593,-0.00003169],
           [450160.398036,-6962890.5431,7.4722,0.007702,-0.00005939]]).dot(np.array([1,T_TDB,T_TDB**2,T_TDB**3,T_TDB**4]));

        for i in range(len(correction3)):
            PHI = sum(correction3[i,0:5]*phi);
            delta_posai=delta_posai+sum(correction3[i,5:7]*np.array([1,T_TDB]))*np.sin(de_ra(PHI/3600.0))+correction3[i,7]*np.cos(de_ra(PHI/3600.0));
            delta_epsilon=delta_epsilon+sum(correction3[i,8:10]*np.array([1,T_TDB]))*np.cos(de_ra(PHI/3600.0))+correction3[i,10]*np.sin(de_ra(PHI/3600.0));
        delta_posai=delta_posai/1e7
        delta_epsilon=delta_epsilon/1e7
        epsilon_bar=84381.406-46.836769*T_TDB-0.0001831*T_TDB**2+0.00200340*T_TDB**3-0.000000576*T_TDB**4-0.0000000434*T_TDB**5;
        # 平黄赤交角

        MOD2TOD = M1(-(epsilon_bar + delta_epsilon)/3600.0).dot(M3(-delta_posai/3600.0).dot(M1(epsilon_bar/3600.0)));
        return MOD2TOD
        #瞬时平天球坐标系to瞬时真天球坐标系

    def TOD_ET(T_TDB,T_TT,DU):

        delta_posai=-0.135e-3*1e7;
        #delta_epsilon = 0.388e-3*1e7;
        phi = np.array([[485868.249036,1717915923.2178,31.8792,0.051635,-0.00024470],
           [1287104.79305,129596581.0481,-0.5532,0.000136,-0.00001149],
           [335779.526232,1739527262.8478,-12.7512,-0.001037,0.00000417],
           [1072260.70369,1602961601.2090,-6.3706,0.006593,-0.00003169],
           [450160.398036,-6962890.5431,7.4722,0.007702,-0.00005939]]).dot(np.array([1,T_TDB,T_TDB**2,T_TDB**3,T_TDB**4]));
        for i in range(len(correction3)):
            PHI=sum(correction3[i,0:5]*phi);
            delta_posai= delta_posai+sum(correction3[i,5:7]*[1,T_TDB])*np.sin(de_ra(PHI/3600.0))+correction3[i,7]*np.cos(de_ra(PHI/3600.0));
            #delta_epsilon = delta_epsilon + sum(correction3[i, 8:10] * np.array([1, T_TDB])) * np.cos(de_ra(PHI / 3600.0)) + correction3[i, 10] * np.sin(de_ra(PHI / 3600.0));
        delta_posai = delta_posai / 1e7
        #delta_epsilon = delta_epsilon/1e7;
        epsilon_bar = 84381.406-46.836769*T_TDB-0.0001831*T_TDB**2+0.00200340*T_TDB**3-0.000000576*T_TDB**4-0.0000000434*T_TDB**5;
        #平黄赤交角

        #-------------------------------------------------------------------------------------------
        theta = 0.7790572732640 + 1.00273781191135448*DU;#地球转过的圈数/圈
        S_bar = 86400*theta+(0.014506+4612.156534*T_TDB+1.3915817*T_TDB**2-0.00000044*T_TDB**3-0.000029956*T_TDB**4-0.0000000368*T_TDB**5)/15;
        #格林尼治平恒星时/秒
        epsilon_gama = delta_posai*np.cos(de_ra(epsilon_bar/3600)) + 0.00264096*np.sin(de_ra(phi[4]/3600))\
                       + 0.00006352*np.sin(de_ra(2*phi[4]/3600))\
                       + 0.00001175*np.sin(de_ra((2*phi[2]-2*phi[3]+3*phi[4])/3600))\
                       + 0.00001121*np.sin(de_ra((2*phi[2]-2*phi[3]+phi[4])/3600))\
                       - 0.00000455*np.sin(de_ra((2*phi[2]-2*phi[3]+2*phi[4])/3600))\
                       + 0.00000202*np.sin(de_ra((2*phi[2]+3*phi[4])/3600)) \
                       + 0.00000198*np.sin(de_ra((2*phi[2]+phi[4])/3600))\
                       - 0.00000172*np.sin(de_ra(3*phi[4]/3600))\
                       - 0.00000087*T_TT*np.sin(de_ra(phi[4]/3600))
        S = S_bar + epsilon_gama/15;
        #格林尼治真恒星时/秒
        TOD2ET = M3(360*S/86400);
        return TOD2ET
        #瞬时真天球坐标系to瞬时地球坐标系
    def ET_CTS(xp,yp):
        ET2CTS=M2(-xp).dot(M1(-yp));
        #瞬时地球坐标系to协议地球坐标系
        return ET2CTS
    #------------------------------------------------------------------------
    #转换矩阵
    CIS2MOD = CIS_MOD(T_TDB); #协议天球坐标系to瞬时平天球坐标系
    MOD2TOD = MOD_TOD(T_TDB); #瞬时平天球坐标系到瞬时真天球坐标系
    TOD2ET = TOD_ET(T_TDB, T_TT, DU); #瞬时真天球坐标系到瞬时地球坐标系
    ET2CTS = ET_CTS(xp, yp); #瞬时地球坐标系到协议地球坐标系

    CIS2CTS=ET2CTS.dot(TOD2ET.dot(MOD2TOD.dot(CIS2MOD)))
    CIS2CTS=np.matrix(CIS2CTS)
    CTS2CIS=CIS2CTS.I
    r_CTS=np.mat(r_CTS)
    r_CTS=r_CTS.I
    #print CIS2CTS

    return (CTS2CIS*r_CTS).I
#print r_CIS([-242691.0,-2101770.5,-6554619.5],2457061.09)