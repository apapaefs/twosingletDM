/* Native regression for the restored historical two-Higgs effective vertices.
 * Copy into a fresh TRSM project and build with
 *   make main=check_restored_operators.c
 * Run with a card and an empty writable TRSM_RUNTIME_DIR.
 * The isolated auxiliary contribution below is a structural regression,
 * not a physical cross section or a full multi-Higgs loop validation.
 */
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

static void require(int condition, const char *message)
{
  if(!condition) { fprintf(stderr,"FAIL: %s\n",message); exit(1); }
}

static void close_value(double actual, double expected, const char *label)
{
  require(isfinite(actual) && isfinite(expected), "nonfinite comparison");
  if(fabs(actual-expected)>1e-10*fabs(expected)+1e-30) {
    fprintf(stderr,"FAIL: %s actual=%.17g expected=%.17g\n",label,actual,expected);
    exit(1);
  }
}

static void refresh(void)
{
  char odd[32];
  require(sortOddParticles(odd)==0, "model parameter evaluation");
}

static double contact(numout *code, double m1, double m2)
{
  /* Fixed 2 TeV total energy, four non-collinear massless outgoing momenta.
   * Only the incoming masses change between the three Higgs pairs.
   */
  double e1=(4000000+m1*m1-m2*m2)/4000, e2=2000-e1;
  double p=sqrt(e1*e1-m1*m1);
  REAL momenta[24]={e1,0,0,p, e2,0,0,-p,
                   500,500,0,0, 500,-500,0,0,
                   500,0,500,0, 500,0,-500,0};
  if(code->interface->nout==2) {
    momenta[8]=1000; momenta[9]=1000;
    momenta[12]=1000; momenta[13]=-1000;
  }
  int error=0;
  require(passParameters(code)==0,"pass process parameters");
  double strong=sqrt(4*acos(-1)*alphaQCD(2000));
  double value=code->interface->sqme(1,strong,momenta,NULL,&error);
  require(error==0 && isfinite(value) && value>0,"nonzero contact contribution");
  return value;
}

int main(int argc, char **argv)
{
  const char *runtime=getenv("TRSM_RUNTIME_DIR");
  require(argc==2 && runtime && *runtime,"provide a card and TRSM_RUNTIME_DIR");
  size_t length=strlen(runtime)+80;
  free(libDir); libDir=malloc(length);
  free(compDir); compDir=malloc(length);
  require(libDir && compDir,"allocate workspace paths");
  snprintf(libDir,length,"%s/so_generated",runtime);
  require(mkdir(libDir,0700)==0,"use a fresh, existing runtime directory");
  snprintf(compDir,length,"%s/comp_%d",runtime,getpid());
  require(readVar(argv[1])==0,"read parameter card");
  require(assignValW("SinT",0.3)==0 && assignValW("Mh2",200)==0,
          "set the native regression point");
  require(assignValW("Maux",1)==0,"restored auxiliary mass parameter");
  refresh();

  double c=findValW("CosT"), s=findValW("SinT");
  double v=2*findValW("MW")*findValW("SW")/findValW("EE");
  double projections[3]={c*c,c*s,s*s};
  double a=alphaQCD(findValW("Mh"))/acos(-1);
  double rqcd=sqrt(1+149.0/12*a+68.6482*a*a-212.447*a*a*a);
  double legacy_aa=-cabs(lAAhiggs(findValW("Mh"),"h1"));
  double legacy_gg=-cabs(lGGhiggs(findValW("Mh"),"h1"))*rqcd;
  char *first[3]={"h1","h1","h2"}, *second[3]={"h1","h2","h2"};
  double reference=0;
  for(int pair=0;pair<3;pair++) {
    double m1=findValW(pair==2 ? "Mh2" : "Mh");
    double m2=findValW(pair==0 ? "Mh" : "Mh2");
    char process[80], library[80];
    for(int gluons=0;gluons<2;gluons++) {
      char *boson=gluons ? "G" : "A";
      snprintf(process,sizeof process,"%s,%s->%s,%s",
               first[pair],second[pair],boson,boson);
      snprintf(library,sizeof library,"restored_contact_%d_%d",pair,gluons);
      numout *code=getMEcode(0,1,process,"h1,h2,A,G",NULL,library);
      require(code && code->interface->nin==2 && code->interface->nout==2,
              "compile the restored two-Higgs/two-vector contact");
      double coupling=gluons ? legacy_gg : legacy_aa;
      /* Vertex: -4*C*p_i*p_j/v (k1.k2*g_mu_nu-k1_nu*k2_mu).
       * Sum polarizations, with 8 final colours for gg; no initial averages.
       * CalcHEP sqme includes the 1/2! identical-final-particle factor.
       */
      double expected=(gluons ? 8 : 1)*4*pow(coupling*projections[pair]*4000000/v,2);
      double value=contact(code,m1,m2);
      close_value(value,expected,"historical two-vector contact normalization");
      printf("PASS %s contact_sqme=%.17g expected=%.17g\n",process,value,expected);
    }
    snprintf(process,sizeof process,"%s,%s->G,G,G,G",first[pair],second[pair]);
    snprintf(library,sizeof library,"restored_aux_%d",pair);
    /* Remove physical exchange diagrams to exercise the x1/G2 chain alone. */
    numout *code=getMEcode(0,1,process,"h1,h2,G",NULL,library);
    require(code && code->interface->nin==2 && code->interface->nout==4,
            "compile the restored two-Higgs/four-gluon chain");
    double value=contact(code,m1,m2);
    if(pair==0) reference=value;
    close_value(value,reference*pow(projections[pair]/projections[0],2),
                "Higgs-pair projection factors");
    require(assignValW("Maux",7)==0,"vary auxiliary mass");
    refresh();
    close_value(contact(code,m1,m2),value,"auxiliary mass independence");
    require(assignValW("Maux",1)==0,"reset auxiliary mass");
    refresh();
    printf("PASS %s auxiliary_sqme=%.17g\n",process,value);
  }
  return 0;
}
