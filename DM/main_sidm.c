/*
 * main_sidm.c -- tree-level CalcHEP evaluation of H3 H3 -> H3 H3 (~X,~X->~X,~X)
 * for cross-checking trsm_sidm.py.  It reuses the TRSM micrOMEGAs project,
 * so the Feynman rules and auto widths are exactly those of the DM pipeline.
 *
 * Build inside an installed TRSM project (see DM/setup_micromegas.sh):
 *     cp DM/main_sidm.c <micromegas>/TRSM/ && make -C <micromegas>/TRSM main=main_sidm.c
 *
 * Usage:
 *     TRSM_RUNTIME_DIR=$(mktemp -d) ./main_sidm card.par [--bins N] [--bwrange R] v1 [v2 ...]
 * where v is the relative velocity in units of c (v = 2 beta_cm, as in
 * trsm_sidm.py).  The card is the usual micrOMEGAs parameter card.
 *
 * Machine-readable output:
 *     TRSM_SIDM_inputs {json}            masses, vev, CalcHEP widths
 *     TRSM_SIDM_cs v=.. pcm=.. sigma_pb=.. err=..
 *     TRSM_SIDM_dcs v=.. cmin=.. cmax=.. sigma_pb=..   (with --bins N)
 *
 * CalcHEP uses fixed widths in s-channel propagators only, and only within
 * its Breit-Wigner range; with --bwrange 1e6 this is exactly mode
 * "tree-fixed-width" in trsm_sidm.py.  cs22 integrates over
 * cos(theta) in [cmin, cmax] and includes the identical-particle factor.
 */
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

void trsm_loop_set_widths(double, double);

static void isolate_workspace(void)
{
  const char *root = getenv("TRSM_RUNTIME_DIR");
  if (!root || !*root) return; /* use the project's own work directory */
  size_t n = strlen(root) + 80;
  free(libDir);  libDir = malloc(n);  snprintf(libDir, n, "%s/so_generated", root); mkdir(libDir, 0700);
  free(compDir); compDir = malloc(n); snprintf(compDir, n, "%s/comp_%d", root, getpid());
}

int main(int argc, char **argv)
{
  int err, nbins = 0, first = 2;
  char cdmName[10];

  double bwrange = -1;
  if (argc < 3) {
    fprintf(stderr, "Usage: %s card.par [--bins N] [--bwrange R] v1 [v2 ...]   (v = relative velocity / c)\n", argv[0]);
    return 2;
  }
  while (first + 1 < argc && strncmp(argv[first], "--", 2) == 0) {
    if (strcmp(argv[first], "--bins") == 0) nbins = atoi(argv[first + 1]);
    else if (strcmp(argv[first], "--bwrange") == 0) bwrange = atof(argv[first + 1]);
    else { fprintf(stderr, "unknown option %s\n", argv[first]); return 2; }
    first += 2;
  }

  ForceUG = 0;
  VZdecay = 1; VWdecay = 1;
  err = readVar(argv[1]);
  if (err == -1) { fprintf(stderr, "Can not open %s\n", argv[1]); return 1; }
  if (err > 0)   { fprintf(stderr, "Wrong card contents at line %d\n", err); return 1; }

  err = sortOddParticles(cdmName);
  if (err) { fprintf(stderr, "Can't calculate %s\n", cdmName); return 1; }

  isolate_workspace();
  double w1 = pWidth("h1", NULL), w2 = pWidth("h2", NULL);
  trsm_loop_set_widths(w1, w2);
  double mx = findValW("MX");
  double vev = 2 * findValW("MW") * findValW("SW") / findValW("EE");
  printf("TRSM_SIDM_inputs {\"Mh\":%.17g,\"Mh2\":%.17g,\"MX\":%.17g,\"vev\":%.17g,"
         "\"SinT\":%.17g,\"CosT\":%.17g,\"vevs\":%.17g,\"LX\":%.17g,\"LHX\":%.17g,\"LSX\":%.17g,"
         "\"width_h1\":%.17g,\"width_h2\":%.17g}\n",
         findValW("Mh"), findValW("Mh2"), mx, vev, findValW("SinT"), findValW("CosT"),
         findValW("vevs"), findValW("LX"), findValW("LHX"), findValW("LSX"), w1, w2);

  numout *cc = newProcess("~X,~X->~X,~X");
  if (!cc) { fprintf(stderr, "newProcess(~X,~X->~X,~X) failed\n"); return 1; }
  passParameters(cc);
  /* CalcHEP keeps the width of an s-channel propagator only while
   * |s - M^2| < ~BWrange*M*Gamma (sqme_aux.c:prepDen) and drops it further
   * away.  --bwrange 1e6 keeps fixed widths everywhere. */
  if (bwrange >= 0) *(cc->interface->BWrange) = bwrange;
  printf("TRSM_SIDM_bwrange %.17g\n", *(cc->interface->BWrange));

  for (int i = first; i < argc; i++) {
    double v = atof(argv[i]);
    if (!(v > 0 && v < 2)) { fprintf(stderr, "skip invalid velocity %s\n", argv[i]); continue; }
    double beta = v / 2, pcm = mx * beta / sqrt(1 - beta * beta);
    int e = 0;
    double sig = cs22(cc, 1, pcm, -1., 1., &e);
    printf("TRSM_SIDM_cs v=%.17g pcm=%.17g sigma_pb=%.17g err=%d\n", v, pcm, sig, e);
    for (int b = 0; b < nbins; b++) {
      double c1 = -1. + 2. * b / nbins, c2 = -1. + 2. * (b + 1) / nbins;
      e = 0;
      double part = cs22(cc, 1, pcm, c1, c2, &e);
      printf("TRSM_SIDM_dcs v=%.17g cmin=%.17g cmax=%.17g sigma_pb=%.17g err=%d\n", v, c1, c2, part, e);
    }
  }
  return 0;
}
