/* Single-scalar SM loops at one loop, GF electroweak input scheme.
 * q-dependent MS effective charm/bottom/strange masses, pole top and W.
 * This helper has no heavy-top QCD K factor or multi-Higgs-loop approximation.
 * The separate historical two-Higgs contacts are defined in the model tables.
 * The on-shell absolute coefficient is sufficient for a partial width;
 * annihilation keeps the complex coefficient and coherent propagators.
 */
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"
/* Exported by SLHAplus in both versions; omitted from the 7.1.4 header. */
extern double MqEff(double mass2GeV, double q);

static double widths[2];
static unsigned long relic_calls=0, indirect_calls=0;
static int loop_stage=0;

static double complex triangle(double tau)
{
  if(tau <= 1) { double a=asin(sqrt(tau)); return a*a; }
  double a=2*acosh(sqrt(tau));
  return -0.25*(a-I*M_PI)*(a-I*M_PI);
}
static double complex fermion(double tau)
{
  if(tau < 1e-4) return 4.0/3+14.0*tau/45+8*tau*tau/63;
  return 2*(tau+(tau-1)*triangle(tau))/(tau*tau);
}
static double complex vector(double tau)
{
  if(tau < 1e-4) return -7-22.0*tau/15-76*tau*tau/105;
  return -(2*tau*tau+3*tau+3*(2*tau-1)*triangle(tau))/(tau*tau);
}

double complex trsm_loop_coefficient(double q, int channel)
{
  double vev=2*findValW("MW")*findValW("SW")/findValW("EE");
  double masses[6]={findValW("Mu"),findValW("Md"),MqEff(0.096,q),McEff(q),MbEff(q),findValW("Mtp")};
  double charges[6]={2.0/3,-1.0/3,-1.0/3,2.0/3,-1.0/3,2.0/3};
  double complex sum=0;
  for(int i=0;i<6;i++) if(masses[i]>0)
    sum += (channel==22 ? 3*charges[i]*charges[i] : 1)*fermion(q*q/(4*masses[i]*masses[i]));
  if(channel==21) return -alphaQCD(q)*sum/(16*M_PI*vev);
  double leptons[3]={0.000510998928,findValW("Mm"),findValW("Ml")};
  for(int i=0;i<3;i++) sum+=fermion(q*q/(4*leptons[i]*leptons[i]));
  double mw=findValW("MW"), ee=findValW("EE");
  sum+=vector(q*q/(4*mw*mw));
  return -ee*ee*sum/(32*M_PI*M_PI*vev);
}

double trsm_loop_abs(double q, double channel)
{ return cabs(trsm_loop_coefficient(q,(int)channel)); }

void trsm_loop_set_widths(double h1, double h2) { widths[0]=h1; widths[1]=h2; }
void trsm_loop_set_stage(int stage) { loop_stage=stage; }
void trsm_loop_report(void)
{ printf("TRSM_loop_hook_v2 relic_calls=%lu indirect_calls=%lu\n",relic_calls,indirect_calls); }

void trsm_loop_probe(double q)
{
  double complex aa=trsm_loop_coefficient(q,22),gg=trsm_loop_coefficient(q,21);
  printf("TRSM_loop_parameters {\"q\":%.17g,\"alpha_s\":%.17g,\"quarks\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],\"leptons\":[%.17g,%.17g,%.17g],\"aa_re\":%.17g,\"aa_im\":%.17g,\"gg_re\":%.17g,\"gg_im\":%.17g}\n",q,alphaQCD(q),findValW("Mu"),findValW("Md"),MqEff(.096,q),McEff(q),MbEff(q),findValW("Mtp"),.000510998928,findValW("Mm"),findValW("Ml"),creal(aa),cimag(aa),creal(gg),cimag(gg));
}

/* Integrated cross section in GeV^-2, including identical final particles.
 * The common SM loop coefficient is multiplied by each doublet projection
 * exactly once. The hXX vertex is twice the potential coefficient K_i33.
 */
double trsm_loop_sigma(double pcm, int channel)
{
  double mx=findValW("MX"), st=findValW("SinT"), ct=findValW("CosT");
  double vev=2*findValW("MW")*findValW("SW")/findValW("EE");
  double q=2*sqrt(mx*mx+pcm*pcm), s=q*q;
  double g1=findValW("LHX")*vev*ct-findValW("LSX")*findValW("vevs")*st;
  double g2=findValW("LHX")*vev*st+findValW("LSX")*findValW("vevs")*ct;
  double m1=findValW("Mh"),m2=findValW("Mh2");
  double complex a=trsm_loop_coefficient(q,channel)*
    (g1*ct/(s-m1*m1+I*m1*widths[0])+g2*st/(s-m2*m2+I*m2*widths[1]));
  double beta=2*pcm/q;
  if(beta<=0) return 0; /* callers use the finite sigma*v limit, pcm>0 */
  return (channel==21 ? 8 : 1)*s*pow(cabs(a),2)/(4*M_PI*beta);
}

void improveCrossSection(long n1,long n2,long n3,long n4,double pcm
#ifdef TRSM_MO7
                         ,double temperature
#endif
                         ,double *result)
{
#ifdef TRSM_MO7
  (void)temperature;
#endif
  if(n1!=pNum("~X") || n2!=pNum("~X") || n3!=n4 || (n3!=21 && n3!=22)) return;
  *result=trsm_loop_sigma(pcm,(int)n3);
  if(loop_stage==1) relic_calls++;
  if(loop_stage==2) indirect_calls++;
}
