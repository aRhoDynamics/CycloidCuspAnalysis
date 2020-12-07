import numpy as np

def stressththe(constA, constB, e, oblq, phase, colat, lon, he, le):
    
    beta20ththe = 0.75*(3.*he - 10.*le)*np.cos(2.*colat) + 0.75*(he - 2.*le)
    beta21ththe = 1.5*(3.*he - 10.*le)*np.sin(2.*colat)
    beta22ththe = -1.5*(3.*he - 10.*le)*np.cos(2.*colat) + 4.5*(he - 2.*le)
    
    ththe = constB*(-6.*e*beta20ththe*np.cos(constA) + e*beta22ththe*(4.*np.sin(2.*lon)*np.sin(constA)+3.*np.cos(2.*lon)*np.cos(constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21ththe*np.cos(lon)*np.sin(phase + constA))
    
    return ththe

def stressphphe(constA, constB, e, oblq, phase, colat, lon, he, le):
    
    beta20phphe = 0.75*(3.*he - 8.*le)*np.cos(2.*colat) + 0.75*(he - 4.*le)
    beta21phphe = 1.5*(3.*he - 8.*le)*np.sin(2.*colat)
    beta22phphe = -1.5*(3.*he - 8.*le)*np.cos(2.*colat) + 4.5*(he - 4.*le)
    
    phphe = constB*(-6.*e*beta20phphe*np.cos(constA) + e*beta22phphe*(4.*np.sin(2.*lon)*np.sin(constA)+3.*np.cos(2.*lon)*np.cos(constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21phphe*np.cos(lon)*np.sin(phase + constA))
    
    return phphe

def stressthphe(constA, constB, e, oblq, phase, colat, lon, he, le):
    
    beta21thphe = 3.*le*np.sin(colat)
    beta22thphe = 3.*le*np.cos(colat)
    
    thphe = constB*(2.*e*beta22thphe*(4.*np.cos(2.*lon)*np.sin(constA)-3.*np.sin(2.*lon)*np.cos(constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21thphe*np.sin(lon)*np.sin(phase + constA))
    
    return thphe

def modeththv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):
    
    beta20ththv = 0.75*(3.*hv - 10.*lv)*np.cos(2.*colat) + 0.75*(hv - 2.*lv)
    beta21ththv = 1.5*(3.*hv - 10.*lv)*np.sin(2.*colat)
    beta22ththv = -1.5*(3.*hv - 10.*lv)*np.cos(2.*colat) + 4.5*(hv - 2.*lv)
    
    ththv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20ththv*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + e*beta22ththv*(4.*np.sin(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda))+3.*np.cos(2.*lon)*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda))) + 4.*np.cos(oblq)*np.sin(oblq)*beta21ththv*np.cos(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))
    
    return ththv

def modephphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):
    
    beta20phphv = 0.75*(3.*hv - 8.*lv)*np.cos(2.*colat) + 0.75*(hv - 4.*lv)
    beta21phphv = 1.5*(3.*hv - 8.*lv)*np.sin(2.*colat)
    beta22phphv = -1.5*(3.*hv - 8.*lv)*np.cos(2.*colat) + 4.5*(hv - 4.*lv)
    
    phphv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20phphv*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + e*beta22phphv*(4.*np.sin(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda))+3.*np.cos(2.*lon)*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda))) + 4.*np.cos(oblq)*np.sin(oblq)*beta21phphv*np.cos(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))
    
    return phphv

def modethphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):
    
    beta21thphv = 3.*lv*np.sin(colat)
    beta22thphv = 3.*lv*np.cos(colat)
    
    thphv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(8.*e*beta22thphv*np.cos(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) - 6.*e*beta22thphv*np.sin(2.*lon)*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21thphv*np.sin(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))
    
    return thphv

def stressththeNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):
    
    alpha22ththe = -1.5*(3.*he_NSR - 10.*le_NSR)*np.cos(2.*colat) + 4.5*(he_NSR - 2.*le_NSR);
    ththeNSR = constB_NSR*alpha22ththe*np.cos(2.*lon + constA_NSR);
    
    return ththeNSR;

def stressphpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):
    
    alpha22phphe = -1.5*(3.*he_NSR - 8.*le_NSR)*np.cos(2.*colat) + 4.5*(he_NSR - 4.*le_NSR)
    phpheNSR = constB_NSR*alpha22phphe*np.cos(2.*lon + constA_NSR)
    
    return phpheNSR

def stressthpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):
    
    alpha22thphe = 3.*le_NSR*np.cos(colat)
    thpheNSR = -2.*constB_NSR*alpha22thphe*np.sin(2.*lon + constA_NSR)
    
    return thpheNSR

def modeththvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):
    
    alpha22ththv = -1.5*(3.*hv_NSR - 10.*lv_NSR)*np.cos(2.*colat) + 4.5*(hv_NSR - 2.*lv_NSR)
    ththvNSR = (1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22ththv*np.cos(2.*lon + constA_NSR - np.arctan(-2.*NSRrate/sj_NSR))
    
    return ththvNSR

def modephphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):
    
    alpha22phphv = -1.5*(3.*hv_NSR - 8.*lv_NSR)*np.cos(2.*colat) + 4.5*(hv_NSR - 4.*lv_NSR)
    phphvNSR = (1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22phphv*np.cos(2.*lon + constA_NSR - np.arctan(-2.*NSRrate/sj_NSR))
    
    return phphvNSR

def modethphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):
    
    alpha22thphv = 3.*lv_NSR*np.cos(colat)
    thphvNSR = -2.*(1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22thphv*np.sin(2.*lon + constA_NSR - np.arctan(-2.*NSRrate/sj_NSR))
    
    return thphvNSR

def getStress(interior, e_in, colat, lon, steps, this_step, oblq, phase, NSRdelta):
    
    import numpy as np
    
    kyear = 1./3.1536E10 # kyear in sec (inverted to ease convertion of s_j from inverse kyears to sec)
    year2sec = 3600.*24.*365.  # conversion for NSR sj values
    
    periodInSec = 306000.0 # sec;
    n = 2.*np.pi/periodInSec # [rad/sec]
    t = (this_step/steps)*periodInSec
    visc = 1E21              # surface viscosity, Pa*s
    ridg = 3.487E9           # surface rigidity, Pa    J-O&V: 3.487E9     Hurford/Rhoden: 3.52E9
    _lambda = ridg/(visc*n)
    sgrav = 1.313            # surface gravity, m/s^2
    radius = 1562.0*1000.    # radius, m               J-O&V: 1562.0      Hurford/Rhoden: 1561.5
    eccCurrent = 0.0094      # eccentricity            J-O&V: 0.0094      Hurford/Rhoden: 0.01
    if (e_in == 1111):
        ecc = eccCurrent
    else:
        ecc = e_in

    if (interior == 1):
        
        numModes = 6
        numModes_NSR = 4
        
        he = 1.15106
        le = 0.308014
        
        # Model 1 values, from J-O&V =VISCO=
        sj1 = kyear*-1.08100E-1 #(*M0*)
        sj2 = kyear*-3.23708E-1 #(*M2*)
        sj3 = kyear*-4.49843E-1 #(*C0*)
        sj4 = kyear*-3.44856   #(*M3*)
        sj5 = kyear*-1.82043E5   #(*T1*)
        sj6 = kyear*-1.02489E6   #(*T2*)
        
        hv1 = 3.69312E-2
        hv2 = 1.38303E-3
        hv3 = 5.17520E-2
        hv4 = 7.16927E-1
        hv5 = 7.19007E-5
        hv6 = 8.49126E-2
        
        lv1 = 1.00431E-2
        lv2 = 2.35498E-4
        lv3 = 1.40299E-2
        lv4 = 1.94722E-1
        lv5 = 6.00150E-3
        lv6 = 2.26341E-2
        
        # These are the s_j values computed by Wade Henning, based on the formulation by J-O&V, for this interior model:
        # Ductile ice viscosity 1e14, Total shell thickness 30km (25km ductile, 5km brittle) 
        # Layer thicknesses: L1=600 km, L2=832 km, L3=100 km, L4=25 km, L5=5 km
        
        he_NSR = 1.90474
        le_NSR = 0.509593
        
        sj1_NSR = -1./(3.37991E6)
        sj2_NSR = -1./(1.1289E6)
        sj3_NSR = -1./(2.00708)
        sj4_NSR = -1./(0.359672)
        
        hv1_NSR = 0.0382668
        hv2_NSR = 0.00204669
        hv3_NSR = 0.000668372
        hv4_NSR = 0.159271
        
        lv1_NSR = 0.0104063
        lv2_NSR = 0.000385665
        lv3_NSR = 0.00968342
        lv4_NSR = 0.0424533

    elif (interior == 2):
    
        numModes = 6
        numModes_NSR = 0
        
        he = 1.151065
        le = 0.3080142
        
        # Model 2 values, from J-O&V =ELASTIC=
        sj1 = kyear*-4.44842E-4
        sj2 = kyear*-1.08230E-1
        sj3 = kyear*-4.49843E-1
        sj4 = kyear*-3.44846
        sj5 = kyear*-1.82134E2
        sj6 = kyear*-1.02492E3
        
        hv1 = 4.31352E-3
        hv2 = 3.43781E-2
        hv3 = 5.15023E-2
        hv4 = 7.17001E-1
        hv5 = 7.10977E-5
        hv6 = 8.47114E-2
        
        lv1 = 1.60251E-1
        lv2 = 9.34250E-3
        lv3 = 1.39915E-2
        lv4 = 1.94817E-1
        lv5 = 5.93813E-3
        lv6 = 2.25805E-2
    
    else:
        print("Undefined interior")
        Quit()

    constA = n*t + np.arctan(_lambda)
    constB = 0.5*((n*n*radius*ridg)/sgrav)*(1./np.sqrt(1. + _lambda*_lambda))

    ththe = stressththe(constA,constB,ecc,oblq,phase,colat,lon,he,le)
    phphe = stressphphe(constA,constB,ecc,oblq,phase,colat,lon,he,le)
    thphe = stressthphe(constA,constB,ecc,oblq,phase,colat,lon,he,le)
    
    if (numModes < 2 or numModes > 6):
        print("Number of modes is disallowed")
        Quit()

    ththv1 = modeththv(n,t,_lambda,sj1,ecc,oblq,phase,colat,lon,hv1,lv1)
    phphv1 = modephphv(n,t,_lambda,sj1,ecc,oblq,phase,colat,lon,hv1,lv1)
    thphv1 = modethphv(n,t,_lambda,sj1,ecc,oblq,phase,colat,lon,hv1,lv1)

    ththv2 = modeththv(n,t,_lambda,sj2,ecc,oblq,phase,colat,lon,hv2,lv2)
    phphv2 = modephphv(n,t,_lambda,sj2,ecc,oblq,phase,colat,lon,hv2,lv2)
    thphv2 = modethphv(n,t,_lambda,sj2,ecc,oblq,phase,colat,lon,hv2,lv2)

    if (numModes > 2):
    
        ththv3 = modeththv(n,t,_lambda,sj3,ecc,oblq,phase,colat,lon,hv3,lv3)
        phphv3 = modephphv(n,t,_lambda,sj3,ecc,oblq,phase,colat,lon,hv3,lv3)
        thphv3 = modethphv(n,t,_lambda,sj3,ecc,oblq,phase,colat,lon,hv3,lv3)
    
    if (numModes > 3):
        
        ththv4 = modeththv(n,t,_lambda,sj4,ecc,oblq,phase,colat,lon,hv4,lv4)
        phphv4 = modephphv(n,t,_lambda,sj4,ecc,oblq,phase,colat,lon,hv4,lv4)
        thphv4 = modethphv(n,t,_lambda,sj4,ecc,oblq,phase,colat,lon,hv4,lv4)
    
    if (numModes > 4):
        
        ththv5 = modeththv(n,t,_lambda,sj5,ecc,oblq,phase,colat,lon,hv5,lv5)
        phphv5 = modephphv(n,t,_lambda,sj5,ecc,oblq,phase,colat,lon,hv5,lv5)
        thphv5 = modethphv(n,t,_lambda,sj5,ecc,oblq,phase,colat,lon,hv5,lv5)
    
    if (numModes > 5):
        
        ththv6 = modeththv(n,t,_lambda,sj6,ecc,oblq,phase,colat,lon,hv6,lv6)
        phphv6 = modephphv(n,t,_lambda,sj6,ecc,oblq,phase,colat,lon,hv6,lv6)
        thphv6 = modethphv(n,t,_lambda,sj6,ecc,oblq,phase,colat,lon,hv6,lv6)

    myStressThTh = ththe + constB*(ththv1 + ththv2 + ththv3 + ththv4 + ththv5 + ththv6)
    myStressPhPh = phphe + constB*(phphv1 + phphv2 + phphv3 + phphv4 + phphv5 + phphv6)
    myStressThPh = thphe + constB*(thphv1 + thphv2 + thphv3 + thphv4 + thphv5 + thphv6) # negative sign comes later

    sum1 = ththv1 + ththv2 + ththv3 + ththv4 + ththv5 + ththv6
    sum2 = phphv1 + phphv2 + phphv3 + phphv4 + phphv5 + phphv6
    sum3 = thphv1 + thphv2 + thphv3 + thphv4 + thphv5 + thphv6
    
    # COMPUTING NSR STRESSES FOR EUROPA
    
    if (NSRdelta != 0):
                
        NSRrate = (ridg/(2.*visc*NSRdelta)) # delta of 0.1 for NSR period of 11419 years; 43 roughly equals 6 Myr period
        NSRperiod = 2.*np.pi/(NSRrate) # testing
        
        constA_NSR = 2.*NSRrate*t + np.arctan(NSRdelta)
        constB_NSR = 0.5*((n*n*radius*ridg)/sgrav)*(1./np.sqrt(1. + NSRdelta*NSRdelta))
                
        ththeNSR = stressththeNSR(constA_NSR,constB_NSR,colat,lon,he_NSR,le_NSR)
        phpheNSR = stressphpheNSR(constA_NSR,constB_NSR,colat,lon,he_NSR,le_NSR)
        thpheNSR = stressthpheNSR(constA_NSR,constB_NSR,colat,lon,he_NSR,le_NSR)
        
        if (numModes_NSR < 2 or numModes_NSR > 4):
            print("Number of NSR modes is disallowed")
            Quit()
    
        ththv1NSR = modeththvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lon,hv1_NSR,lv1_NSR)
        phphv1NSR = modephphvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lon,hv1_NSR,lv1_NSR)
        thphv1NSR = modethphvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lon,hv1_NSR,lv1_NSR)
        
        ththv2NSR = modeththvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lon,hv2_NSR,lv2_NSR)
        phphv2NSR = modephphvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lon,hv2_NSR,lv2_NSR)
        thphv2NSR = modethphvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lon,hv2_NSR,lv2_NSR)
        
        if (numModes_NSR > 2):
            
            ththv3NSR = modeththvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lon,hv3_NSR,lv3_NSR)
            phphv3NSR = modephphvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lon,hv3_NSR,lv3_NSR)
            thphv3NSR = modethphvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lon,hv3_NSR,lv3_NSR)
        
        if (numModes_NSR > 3):
            
            ththv4NSR = modeththvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lon,hv4_NSR,lv4_NSR)
            phphv4NSR = modephphvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lon,hv4_NSR,lv4_NSR)
            thphv4NSR = modethphvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lon,hv4_NSR,lv4_NSR)
        
        # Combining stresses

        myStressThThNSR = ththeNSR + constB_NSR*(ththv1NSR + ththv2NSR + ththv3NSR + ththv4NSR)
        myStressPhPhNSR = phpheNSR + constB_NSR*(phphv1NSR + phphv2NSR + phphv3NSR + phphv4NSR)
        myStressThPhNSR = thpheNSR + constB_NSR*(thphv1NSR + thphv2NSR + thphv3NSR + thphv4NSR)
        
        myStressThThTot = myStressThTh + myStressThThNSR
        myStressPhPhTot = myStressPhPh + myStressPhPhNSR
        myStressThPhTot = myStressThPh + myStressThPhNSR

    else:
        myStressThThTot = myStressThTh
        myStressPhPhTot = myStressPhPh
        myStressThPhTot = myStressThPh
    
    # For writing out (west longitude) only
    latDeg = 90-np.degrees(colat)
    lonDeg = 360. - np.degrees(lon)
    meanMotion = 2.*np.pi*(this_step/steps)
    
    # Converting COMBINED STRESSES to principle stresses.
    
    if (myStressThThTot == myStressPhPhTot):
        zeta = np.pi/2
    else:
        zeta = 0.5*np.arctan((2.*myStressThPhTot)/(myStressThThTot-myStressPhPhTot)) # Amount by which coordinates are being rotated to get principal stresses; changes CCW, orientation of sigTheta

    sigTheta   = myStressThThTot*(np.square(np.cos(zeta)))+myStressPhPhTot*(np.square(np.sin(zeta)))+myStressThPhTot*np.sin(2.*zeta)  # Corresponds to Hurford's sigma 1
    sigPhi = myStressThThTot*(np.square(np.sin(zeta)))+myStressPhPhTot*(np.square(np.cos(zeta)))-myStressThPhTot*np.sin(2.*zeta)      # Corresponds to Hurford's sigma 2

    sigThetaKPA = sigTheta*1E-3
    sigPhiKPA   = sigPhi*1E-3
    
    zetaCW = (2.*np.pi)-zeta   # Changes to CW, still oriented along sigTheta
    
    # stress: The largest (i.e. most tensile) principal stress
    # heading: Direction of crack propagation/formation, perpendicular to direction of largest principal stress (i.e. max tension)
    
    # Determines which stress is most tensile and largest heading.
    
    if (sigThetaKPA < sigPhiKPA):
        stress = sigPhiKPA         # In this case, sigPhi is the max stress, so the heading should be perpendicular to sigPhi. SigTheta is -| to sigPhi by definition, so heading = zetaCW
        heading = zetaCW
    else:
        stress = sigThetaKPA
        heading = zetaCW +(np.pi/2.)  # In this case, sigTheta is the max stress, so the heading should be perpendicular to sigTheta. SigTheta is oriented along zeta, so heading = zetaCW + 90 deg
    
    if (heading >= (2.*np.pi)):       # Making sure azimuth values fall between 0 and 360. 
        heading = heading -(2.*np.pi) # Also finds the two, complementary heading directions (e.g. 45 and 135)
        heading2 = heading +np.pi
    else:
        heading2 = heading - np.pi
    
    if (heading > heading2):       # Determines which of the two heading directions (e.g. 45 or 135) is largest and selects that as the output heading.
        bigHeading = heading
    else:
        bigHeading = heading2
        
    return (stress, bigHeading)
