#include <string.h>
#include <sys/stat.h>
#include <unistd.h>


/*====== Modules ===============
   Keys to switch on
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO
  /* Display information about mass spectrum  */

#define CONSTRAINTS
//#define SMODELS
//#define MONOJET
//#define HIGGSBOUNDS
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS

#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */

#define INDIRECT_DETECTION
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation;
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line
     of sight;
     Calculate galactic propagation of positrons and antiprotons.
  */

//#define RESET_FORMFACTORS
  /* Modify default nucleus form factors,
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */
#define CDM_NUCLEON
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */

// #define CDM_NUCLEUS
     // Calculate  exclusion rate for direct detection experiments Xenon1T, DarkSide50, CRESST, and PICO

//#define NEUTRINO
 /*  Neutrino signal of DM annihilation in Sun and Earth */

// #define DECAYS

//#define CROSS_SECTIONS

/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS
     /* Display  graphical plots on the screen */
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"



void trsm_loop_set_widths(double,double);
void trsm_loop_set_stage(int);
void trsm_loop_report(void);
double trsm_loop_sigma(double,int);
void trsm_loop_probe(double);

static void isolate_loop_workspace(void)
{
  const char *root=getenv("TRSM_RUNTIME_DIR");
  if(!root || !*root) { fprintf(stderr,"TRSM_RUNTIME_DIR is required for isolated v2 execution\n"); exit(2); }
  size_t n=strlen(root)+80;
  free(libDir); libDir=malloc(n); snprintf(libDir,n,"%s/so_generated",root); mkdir(libDir,0700);
  free(compDir); compDir=malloc(n); snprintf(compDir,n,"%s/comp_%d",root,getpid());
}

int main(int argc,char** argv)
{  int err;
   int planckCMBEnabled=0;
   char cdmName[10];
   int spin2, charge3,cdim;

   /* This query must not load a parameter card or modify the CalcHEP workspace. */
   if(argc==2 && strcmp(argv[1],"--capabilities")==0)
   {
     int virtual_enabled=getenv("TRSM_LEGACY_VIRTUAL_OFF")==NULL;
     printf("{\"schema\":\"trsm_driver_capabilities_v1\",\"planck_cmb\":{"
            "\"method\":\"micromegas_planck2018_swave_v1\","
            "\"pann_limit_cm3_s_GeV\":3.2e-28,\"spectrum_key\":7,"
            "\"spectra_flag\":%d,\"vrot_km_s\":%.17g,"
            "\"vz_decay\":%d,\"vw_decay\":%d},\"physics_version\":\"trsm_constraints_v2\",\"loop_hook\":true}\n",SpectraFlag,vRot,virtual_enabled,virtual_enabled);
     return 0;
   }
   if(argc>2)
   {
     if(argc!=3 || strcmp(argv[2],"--planck-cmb")!=0)
     { fprintf(stderr,"Usage: %s parameter-card [--planck-cmb]\n",argv[0]); return 2; }
     planckCMBEnabled=1;
   }

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  VZdecay=1; VWdecay=1;
  if(getenv("TRSM_LEGACY_VIRTUAL_OFF")) { VZdecay=0; VWdecay=0; }

  if(argc==1)
  {
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }

  err=readVar(argv[1]);

  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}



  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  isolate_loop_workspace();
  double trsm_width_h1=pWidth("h1",NULL), trsm_width_h2=pWidth("h2",NULL);
  trsm_loop_set_widths(trsm_width_h1,trsm_width_h2);
  printf("TRSM_inputs_v2 {\"Mh\":%.17g,\"Mh2\":%.17g,\"MX\":%.17g,\"EE\":%.17g,\"SW\":%.17g,\"MW\":%.17g,\"GF\":%.17g,\"Mtp\":%.17g,\"MbMb\":%.17g,\"McMc\":%.17g,\"alphaSMZ\":%.17g,\"width_h1\":%.17g,\"width_h2\":%.17g,\"VWdecay\":%d,\"VZdecay\":%d}\n",
    findValW("Mh"),findValW("Mh2"),findValW("MX"),findValW("EE"),findValW("SW"),findValW("MW"),findValW("GF"),findValW("Mtp"),findValW("MbMb"),findValW("McMc"),findValW("alphaSMZ"),trsm_width_h1,trsm_width_h2,VWdecay,VZdecay);
  if(getenv("TRSM_LOOP_PROBE")) {
    double pcm=atof(getenv("TRSM_LOOP_PROBE"));
    double mx=findValW("MX");
    trsm_loop_probe(2*sqrt(mx*mx+pcm*pcm));
    for(int i=1;i<=2;i++) {
      char name[8]; sprintf(name,"h%d",i); txtList channels=NULL;
      double width=pWidth(name,&channels);
      printf("TRSM_partial_width h%d aa=%.17g gg=%.17g\n",i,width*findBr(channels,"A,A"),width*findBr(channels,"G,G"));
      trsm_loop_probe(findValW(i==1?"Mh":"Mh2"));
    }
    printf("TRSM_loop_probe pcm=%.17g aa=%.17g gg=%.17g\n",pcm,trsm_loop_sigma(pcm,22),trsm_loop_sigma(pcm,21));
    if(getenv("TRSM_CALCHEP_PROBE")) {
      numout *code=newProcess("~X,~X->A,A"); int error=0;
      if(code) { passParameters(code); printf("TRSM_calchep_probe aa_pb=%.17g\n",cs22(code,1,pcm,-1,1,&error)); }
      code=newProcess("~X,~X->G,G"); error=0;
      if(code) { passParameters(code); printf("TRSM_calchep_probe gg_pb=%.17g\n",cs22(code,1,pcm,-1,1,&error)); }
    }
    return 0;
  }
  for(int k=1;k<=Ncdm;k++)
  {
     qNumbers(CDM[k], &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.17E\n",CDM[k],  spin2,McdmN[k]);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.17E pb\n",csLim);
}
#endif

#ifdef SMODELS
{
  int combineSRs=0;
  char* combineAnas=NULL;//"ATLAS-SUSY-2018-05-ewk,ATLAS-SUSY-2019-08,ATLAS-SUSY-2019-09";

  int status=0, smodelsOK=0;
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  double CombRvalue, CombRexpected, CombSmoLsig, CombSmoLmax, CombSmoLSM;

  char analysis[50]={},topology[100]={},smodelsInfo[100];
  char CombAnalyses[200]={};
  int LHCrun=LHC8|LHC13;  //  LHC8 - 8TeV; LHC13 - 13TeV;
//  int LHCrun=LHC13;  //  LHC13 - 13TeV only;

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK)
  { printf("\n highest r-value = %.17E",Rvalue);

    if(Rvalue>0)
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0)
      { printf("\n expected r = %.17E ",Rexpected);
        if(SmoLsig!=INFINITY)
        { printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.17E %.17E",
                  2*(SmoLsig-SmoLmax), 2*(SmoLSM-SmoLmax) );
        }
      }
    }
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded");
    else if(status==-1) printf("\n not not tested by results in SModelS database");
    printf("\n");

    // r-value and likelihoods from analysis cvombination
    if(CombRvalue>0)
    { printf("\n Combination of %s",CombAnalyses);
      printf("\n r-value = %.17E (expected r = %.17E)",CombRvalue, CombRexpected);
      if(CombRvalue>=1) printf("  --> excluded");
      else printf("  --> not excluded");
      printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.17E %.17E \n\n",
                    2*(CombSmoLsig-CombSmoLmax),2*(CombSmoLSM-CombSmoLmax));
    }

  } else system("cat smodels.err"); // problem: see smodels.err
}

#endif


#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);
}
#endif

#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50];
   NH0=hbBlocksMO("HB.in",&NHch);
//    NH0= hbBlocksMDL("HB.in",&NHch);
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.17E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version);
     printf("  Nobservables=%.0f chi^2 = %.17E pval= %.17E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char Lilith_version[50];
   if(LilithMO("Lilith_in.xml"))
   {
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.17E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }
   } else printf("LILITH: there is no Higgs candidate\n");
}
#endif


#ifdef SMODELS
{ int combineSRs=0;
  char* combineAnas=NULL;

  int status=0, smodelsOK=0;
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  double CombRvalue, CombRexpected, CombSmoLsig, CombSmoLmax, CombSmoLSM;

  char analysis[50]={},topology[100]={},smodelsInfo[100];
  char CombAnalyses[200]={};
  int LHCrun=LHC8|LHC13;  //  LHC8 - 8TeV; LHC13 - 13TeV;
//  int LHCrun=LHC13;  //  LHC13 - 13TeV only;

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK)
  { printf("\n highest r-value = %.17E",Rvalue);

    if(Rvalue>0)
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0)
      { printf("\n expected r = %.17E ",Rexpected);
        if(SmoLsig!=INFINITY)
        { printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.17E %.17E",
                  2*(SmoLsig-SmoLmax), 2*(SmoLSM-SmoLmax)); }
      }
    }
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded");
    else if(status==-1) printf("\n not not tested by results in SModelS database");
    printf("\n");

    // r-value and likelihoods from analysis cvombination
    if(CombRvalue>0)
    { printf("\n Combination of %s",CombAnalyses);
      printf("\n r-value = %.17E (expected r = %.17E)",CombRvalue, CombRexpected);
      if(CombRvalue>=1) printf("  --> excluded");
      else printf("  --> not excluded");
      printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.17E %.17E \n\n",
                    -2*(CombSmoLsig-CombSmoLmax),-2*(CombSmoLSM-CombSmoLmax));
    }

  } else system("cat smodels.err"); // problem: see smodels.err
}

#endif


#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega;
  int i,err;
  printf("\n==== Calculation of relic density =====\n");

  if(Ncdm==1)
  {  double Xf;
     trsm_loop_set_stage(1);
     Omega=darkOmega(&Xf,fast,Beps,&err);
     trsm_loop_set_stage(0);
     printf("Xf=%.17g Omega=%.17g darkOmega_error=%d\n",Xf,Omega,err);
     if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  } else
  if(Ncdm==2)
  {
    Omega= darkOmega2(fast,Beps,&err);
    printf("Omega_1h^2=%.17E Omega_2h^2=%.17E err=%d \n", Omega*fracCDM[1], Omega*fracCDM[2],err);
  }else
  {
     Omega=darkOmegaN(fast,Beps,&err);
     printf("Omega=%.17E\n",Omega);
     for(int k=1;k<=Ncdm;k++) printf("   Omega_%d=%.17E\n",k,Omega*fracCDM[k]);
  }
}

#endif

#ifdef FREEZEIN
{
  double TR=1E6;
  double omegaFi;
  toFeebleList(CDM[1]);
  VWdecay=0; VZdecay=0;

  omegaFi=darkOmegaFi(TR,CDM[1],&err);
  printf("omega freeze-in=%.17E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}
#endif



#ifdef INDIRECT_DETECTION
trsm_loop_set_stage(2);
{
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;

printf("\n==== Indirect detection =======\n");

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.

       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections
    */
  if(planckCMBEnabled)
  {
    const char *reason=NULL;
    double cmbRatio=NAN;
    if(err) reason="calcSpectrum_error";
    else if(!isfinite(sigmaV) || sigmaV<0) reason="invalid_annihilation_rate";
    else
    {
      for(i=0;i<NZ;i++)
        if(!isfinite(SpA[i]) || !isfinite(SpE[i])) reason="nonfinite_spectrum";
      if(SpA[0]<=0 || SpE[0]<=0) reason="invalid_spectrum_mass";
      if(!reason)
      {
        cmbRatio=PlanckCMB(sigmaV,SpA,SpE);
        if(!isfinite(cmbRatio) || cmbRatio<0) reason="invalid_cmb_ratio";
      }
    }
    if(reason)
      printf("TRSM_PlanckCMB_v1 {\"status\":\"error\",\"reason\":\"%s\","
             "\"spectrum_error\":%d,\"ratio_raw\":null}\n",reason,err);
    else
      printf("TRSM_PlanckCMB_v1 {\"status\":\"ok\",\"ratio_raw\":%.17g,"
             "\"sigma_v_cm3_s\":%.17g}\n",cmbRatio,sigmaV);
  }
  /* --- Photon continuum diagnostic: not used for Fermi-LAT line exclusion --- */
  {
    double fi = 0.1, dfi = 0.05; /* old diagnostic cone */
    gammaFluxTab(fi, dfi, sigmaV, SpA, FluxA);

    printf("Photon continuum flux diagnostic for angle of sight f=%.2f[rad]\n"
          "and spherical region described by cone with angle %.2f[rad]\n",
          fi, 2*dfi);

  #ifdef SHOWPLOTS
    sprintf(txt, "Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",
            fi, 2*dfi);
    displayPlot(txt, "E[GeV]", Emin, Mcdm, 0, 1, "", 0, SpectdNdE, FluxA);
  #endif

    printf("Photon continuum differential flux = %.17E[cm^-2 s^-1 GeV^-1] "
          "for E=%.1f[GeV]\n",
          SpectdNdE(Etest, FluxA), Etest);
    printf("This differential continuum value is not used for Fermi-LAT line limits.\n");
  }

  /* --- Fermi-LAT gamma-line observable --- */
  {
    int ch;
    int foundLine = 0;
    int invalidLine = 0;
    double totalLineFluxR16 = 0.0;

    const double deg = M_PI/180.0;
    const double rGC = 16.0*deg;       /* Fermi R16 radius */
    const double maskB = 5.0*deg;      /* Galactic-plane mask |b| < 5 deg */
    const double maskL = 6.0*deg;      /* mask applies only for |l| > 6 deg */
    const double step = 0.5*deg;       /* numerical patch size for ROI sum */

    printf("\nFermi-LAT gamma-line flux estimate using R16 ROI approximation\n");
    printf("ROI: sqrt(l^2+b^2)<16 deg, excluding |b|<5 deg and |l|>6 deg\n");

    if (!vSigmaCh) {
      printf("FermiLAT_line_flux_R16 could not be calculated properly: "
            "vSigmaCh is not available after calcSpectrum.\n");
    } else {
      for (ch = 0; vSigmaCh[ch].weight > 0; ++ch) {
        char *p3 = vSigmaCh[ch].prtcl[2];
        char *p4 = vSigmaCh[ch].prtcl[3];
        char *p5 = vSigmaCh[ch].prtcl[4];

        int nGamma = 0;
        char *xParticle = NULL;
        double sigmaVChannel;
        double sigmaVPhotonWeighted;
        double mX = 0.0;
        double eGamma;
        double lineFlux = 0.0;

        double l, b;

        /* Fermi-LAT line limits are for monochromatic 2-body lines.
          Skip 3-body gamma radiation such as XX -> A W+ W-. */
        if (p5 != NULL) continue;
        if (p3 == NULL || p4 == NULL) continue;

        if (strcmp(p3, "A") == 0) nGamma++;
        else xParticle = p3;

        if (strcmp(p4, "A") == 0) nGamma++;
        else xParticle = p4;

        if (nGamma == 0) continue;

        foundLine = 1;

        sigmaVChannel = sigmaV * vSigmaCh[ch].weight;
        sigmaVPhotonWeighted = nGamma * sigmaVChannel;

        if (nGamma == 1) {
          if (xParticle == NULL) {
            invalidLine = 1;
            printf("Skipping photon-line channel %s %s: missing recoil particle.\n",
                  p3, p4);
            continue;
          }

          mX = pMass(xParticle);
          if (mX <= 0.0 || 2.0*Mcdm <= mX) {
            invalidLine = 1;
            printf("Skipping photon-line channel %s %s: invalid recoil mass "
                  "mX=%.17E GeV for mDM=%.17E GeV.\n",
                  p3, p4, mX, Mcdm);
            continue;
          }

          eGamma = Mcdm * (1.0 - mX*mX/(4.0*Mcdm*Mcdm));
        } else if (nGamma == 2) {
          eGamma = Mcdm;
        } else {
          invalidLine = 1;
          printf("Skipping photon-line channel %s %s: unexpected photon count %d.\n",
                p3, p4, nGamma);
          continue;
        }

        if (eGamma < 0.214 || eGamma > 462) {
          printf("Photon-line channel %s %s has E_gamma=%.17E GeV, outside "
                "Fermi-LAT table range 0.214-462 GeV. Flux not printed for limits.\n",
                p3, p4, eGamma);
          continue;
        }

        /* Sum gammaFluxGC over small Galactic-coordinate rectangles to approximate
          the Fermi R16 ROI mask. gammaFluxGC returns integrated flux [cm^-2 s^-1]
          when passed an integrated line cross section [cm^3 s^-1]. */
        for (l = -rGC; l < rGC; l += step) {
          for (b = -rGC; b < rGC; b += step) {
            double lc = l + 0.5*step;
            double bc = b + 0.5*step;

            if (sqrt(lc*lc + bc*bc) > rGC) continue;
            if (fabs(bc) < maskB && fabs(lc) > maskL) continue;

            lineFlux += gammaFluxGC(l, b, step, step, sigmaVPhotonWeighted);
          }
        }

        totalLineFluxR16 += lineFlux;

        printf("FermiLAT_line_channel %s %s: E_gamma=%.17E[GeV], "
              "sigmaV=%.17E[cm^3 s^-1], N_gamma*sigmaV=%.17E[cm^3 s^-1], "
              "Phi_R16=%.17E[cm^-2 s^-1]\n",
              p3, p4, eGamma, sigmaVChannel, sigmaVPhotonWeighted, lineFlux);
      }

      if (!foundLine) {
        printf("FermiLAT_line_flux_R16 could not be calculated properly: "
              "no 2-body photon-line channels were found in vSigmaCh. "
              "The printed continuum photon flux must not be compared to "
              "Fermi-LAT line limits.\n");
      } else if (totalLineFluxR16 > 0.0) {
        printf("FermiLAT_total_line_flux_R16 = %.17E[cm^-2 s^-1]\n",
              totalLineFluxR16);
      } else {
        printf("FermiLAT_line_flux_R16 could not be calculated properly: "
              "line channels were found, but none were inside the usable "
              "Fermi-LAT 0.214-462 GeV energy range or passed validity checks.\n");
      }

      if (invalidLine) {
        printf("Warning: at least one photon-line channel failed validity checks.\n");
      }
    }
  }



// //  if(SpA)
//   {
//      double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */

//      gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);
//      printf("Photon flux  for angle of sight f=%.2f[rad]\n"
//      "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
// #ifdef SHOWPLOTS
//      sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
//      displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
// #endif
//      printf("Photon flux = %.17E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);
//   }

// //  if(SpE)
//   {
//     posiFluxTab(Emin, sigmaV, SpE,  FluxE);
// #ifdef SHOWPLOTS
//     displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
// #endif
//     printf("Positron flux  =  %.17E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
//     SpectdNdE(Etest, FluxE),  Etest);
//   }

// //  if(SpP)
//   {
//     pbarFluxTab(Emin, sigmaV, SpP,  FluxP  );
// #ifdef SHOWPLOTS
//      displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
// #endif
//     printf("Antiproton flux  =  %.17E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
//     SpectdNdE(Etest, FluxP),  Etest);
//   }
}
#endif

#ifdef RESET_FORMFACTORS
{
/*
   The user has approach to form factors  which specifies quark contents
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma";
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");

  printf("protonFF (default) d %.17E, u %.17E, s %.17E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(default) d %.17E, u %.17E, s %.17E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
//                    To restore default form factors of  version 2  call
     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.17E, u %.17E, s %.17E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(new)     d %.17E, u %.17E, s %.17E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call
//  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;
  double csSIp1,csSIn1,csSDp1,csSDn1, csSIp1_,csSIn1_,csSDp1_,csSDn1_;
  double csSIp2,csSIn2,csSDp2,csSDn2, csSIp2_,csSIn2_,csSDp2_,csSDn2_;
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");

  for(int k=1;k<=Ncdm;k++)
  {
    nucleonAmplitudes(CDM[k], pA0,pA5,nA0,nA5);
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM[k],antiParticle(CDM[k]));
    printf("proton:  SI  %.17E [%.17E]  SD  %.17E [%.17E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.17E [%.17E]  SD  %.17E [%.17E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] );

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*McdmN[k]/(Nmass+ McdmN[k]),2.);
    csSIp1=  SCcoeff*pA0[0]*pA0[0];  csSIp1_=  SCcoeff*pA0[1]*pA0[1];
    csSDp1=3*SCcoeff*pA5[0]*pA5[0];  csSDp1_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn1=  SCcoeff*nA0[0]*nA0[0];  csSIn1_=  SCcoeff*nA0[1]*nA0[1];
    csSDn1=3*SCcoeff*nA5[0]*nA5[0];  csSDn1_=3*SCcoeff*nA5[1]*nA5[1];

    printf("%s[%s]-nucleon cross sections[pb]:\n",CDM[k],antiParticle(CDM[k]));
    printf(" proton  SI %.17E [%.17E] SD %.17E [%.17E]\n", csSIp1,csSIp1_,csSDp1,csSDp1_);
    printf(" neutron SI %.17E [%.17E] SD %.17E [%.17E]\n", csSIn1,csSIn1_,csSDn1,csSDn1_);
  }
}
#endif

#ifdef CDM_NUCLEUS
{ char* expName;
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
       if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval));
  else printf("Not excluded by DD experiments  at 90%% level \n");
}
#endif

#ifdef NEUTRINO
if(!CDM[1] || !CDM[2])
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;

 printf("\n===============Neutrino Telescope=======  for  ");
 if(forSun) printf("Sun\n"); else printf("Earth\n");

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif
{
    printf(" E>%.1E GeV neutrino flux       %.17E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.17E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));
}

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);
#endif
    printf(" E>%.1E GeV Upward muon flux    %.17E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu);
#endif
  printf(" E>%.1E GeV Contained muon flux %.17E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}
#endif


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width;
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  }

}
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;

  next=nextOdd(1,&nextM);
  if(next && nextM<1000)
  {
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.17E GeV\n",2*Pcm);

     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.17E[pb]\n",next,next_, cs);
  }
}

#endif

#ifdef CLEAN
trsm_loop_report();
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");
#endif



  killPlots();
  return 0;
}
