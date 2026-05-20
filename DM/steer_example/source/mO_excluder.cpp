#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <cerrno>
#include <vector>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>

using namespace std;

namespace {

bool ensure_directory(const string& path) {
  if (mkdir(path.c_str(), 0777) == 0) {
    return true;
  }

  return errno == EEXIST;
}

string output_path(const string& directory, const string& filename) {
  if (directory.empty() || directory == ".") {
    return filename;
  }
  return directory + "/" + filename;
}

void append_line(const string& path, const string& line) {
  ofstream stream(path.c_str(), ios::app);
  stream << line << endl;
}

string shell_quote(const string& text) {
  string quoted = "'";
  for (string::const_iterator it = text.begin(); it != text.end(); ++it) {
    if (*it == '\'') {
      quoted += "'\\''";
    } else {
      quoted += *it;
    }
  }
  quoted += "'";
  return quoted;
}

}  // namespace


/*
double LUXexcl(double mH){
double excl;
// for the fitted function see LUX_fit
excl = 1.86934E-10 + 1.1236E-11*mH + (-1.32312E-16)*mH*mH + (1.85167E8)/((-9.35439E15)*mH + 1.12151E15*mH*mH);
return excl;}

*/

// double LUXexcl(double mH){                                                                                                                                                 double excl;

//   if ((mH < 17.) || mH > 200.)
//     excl=-0.138541 + 0.00343701*mH - 1.26198e-6*mH*mH + 3.44392e+10/(-1.49313e+9*mH + 2.32828e+8*mH*mH);
//   else if (mH < 33.)
//     excl=-0.0209648 - 0.00479147*mH + 0.0000388817*mH*mH + 2.75483/(0.248816*mH - 0.00173454*mH*mH);
//   else excl=0.295076 - 0.00150593*mH + 0.0000109119*mH*mH + 0.999075/(1.00003*mH + 1.0009*mH*mH);
//   excl=excl*1.e-9;
//   return excl;}

double DirDetexcl(double x){

  double result;
  if (x<10){
    result= (exp(-23.2949 + 2.2493* pow((-3.60228 + log(x)),2)));
    result=result*(1-0.465188 + 1.64349*pow((-2.3605 + log(x)),2));}
  // // else if (x<11) result=(8.75459*1.e-9 - 6.96272*1.e-10* x);
  // //  else if (x<18) {result=exp(-23.2938 + 0.0135204* pow((-24.8742 + x),2));
  // //  result=result*(1-0.0826054 + 0.0110134*pow((-13.3867 + x),2));}
  // // else if (x<20)result= 5.0606*1.e-10 - 1.91361*1.e-11* x;
  // else if (x<20) result=(9.301626524677054e-8*exp(33.018821405901335/log(x))*pow(x,1.5351327504595669))*1.e-10;
  // // else if (x<30) result=5.103100000000005e-10 - 3.674166666666671e-11*x + 7.027333333333344e-13*pow(x,2);
  // else if (x<30) result=1.744393487921861e-28*exp(68.14955833703276/log(x))*pow(x,5.8655281320921215);

  // //  else if (x<30) result=(7.64124*1.e-11 + 1.34496*1.e-13*pow ((-36.6817 + x),2));

  else if (x<40) result=2.45 * 1.e-12 * pow(10, (1.65 / pow(1 - log10(40),2)) * pow(log10(x) - log10(40),2)  ); // update LZ2025
  // else if (x<40) result= (2.771+0.0427*x)*1.e-11;
  // else if (x<40) result= (3.371-0.0227*x)*1.e-12; // update LZ2025
  // //  else if (x<65) result=3.050728*1.e-13*x-3.13092*1.e-12; // new lux zeplin for this region only
  // else if (x< 60) result= 4.84775*1.e-11 + 7.28365*1.e-13* x;
  else if (x< 60) result= 2.32775*1.e-12 + 2.28365*1.e-15* x; // update LZ2025
  else if (x< 80) result= 1.68775*1.e-12 + 1.28365*1.e-14* x; // update LZ2025

  // // else if (x<100) result = exp(0.671326*log(x)-25.8196);
  // // else result=1.49167*1.e-12*exp(0.977674*log(x));
  else if (x<100) result = exp(0.671326*log(x)-29.5696); // update LZ2025
  else result=3.49167*1.e-14*exp(0.977674*log(x)); // update LZ2025

  //  if (( x> 400) && (result > 8.243e-10)) result=8.243e-10; //very rough cutoff for new xenon1t

  // if ((x> 40) && (x<100)) result=(0.0774*x+1.382)*1.e-11;

  // if ((x> 50) && (x<70)) result=(0.0767*x+1.423)*1.e-11;
  // else if ((x > 100) && (x <200)) result = (0.0805*x+1.235)*1.e-11;
  // else if ((x>200) && (x<400)) result=8.135e-13*x+8.8e-12;
  // else if (x>400) result=8.168e-13*x+7.5e-12; // very rough approx for new xenon1t

  // if ((x>40) && (x<65)) result=3.050728*1.e-13*x-3.13092*1.e-12; // new lux zeplin for this region only
  // if ((x>= 65) && (x<110)) result=2.80493*1.e-13*x-1.53327*1.e-12; // same
  // if ((x>=78) && (x < 219)) result=(0.0281889*x-0.154847)*1.e-11; // same

  return result;


  //  double result = exp(0.671326*log(mH)-25.8196);
  // return result;
}



/*


double DirDetexcl(double x){

  double result;
  if (x<10){
    result= (exp(-23.2949 + 2.2493* pow((-3.60228 + log(x)),2)));
    result=result*(1-0.465188 + 1.64349*pow((-2.3605 + log(x)),2));}
  else if (x<11) result=(8.75459*1.e-9 - 6.96272*1.e-10* x);
  else if (x<18) {result=exp(-23.2938 + 0.0135204* pow((-24.8742 + x),2));
    result=result*(1-0.0826054 + 0.0110134*pow((-13.3867 + x),2));}
  else if (x<20)result= 5.0606*1.e-10 - 1.91361*1.e-11* x;
  else if (x<38) result=(7.64124*1.e-11 + 1.34496*1.e-13*pow ((-36.6817 + x),2));
  else if (x< 60) result= 4.84775*1.e-11 + 7.28365*1.e-13* x;

  else if (x<100) result = exp(0.671326*log(x)-25.8196);
  else result=1.49167*1.e-12*exp(0.977674*log(x));

  //  if (( x> 400) && (result > 8.243e-10)) result=8.243e-10; //very rough cutoff for new xenon1t

  if ((x> 40) && (x<100)) result=(0.0774*x+1.382)*1.e-11;

  if ((x> 50) && (x<70)) result=(0.0767*x+1.423)*1.e-11;
  else if ((x > 100) && (x <200)) result = (0.0805*x+1.235)*1.e-11;
  else if ((x>200) && (x<400)) result=8.135e-13*x+8.8e-12;
  else if (x>400) result=8.168e-13*x+7.5e-12; // very rough approx for new xenon1t
  return result;


  //  double result = exp(0.671326*log(mH)-25.8196);
  // return result;
}


*/


void write_dirdet_limit_plot(const string& outputDir) {
  if (!ensure_directory(outputDir)) {
    cerr << "Could not create output directory " << outputDir << endl;
    return;
  }

  const string dataPath = output_path(outputDir, "a_dirdet_limit_curve.dat");
  const string scriptPath = output_path(outputDir, "a_dirdet_limit_curve.gnuplot");
  const string plotPath = output_path(outputDir, "a_dirdet_limit_curve.png");

  ofstream data(dataPath.c_str());
  data << "# mDM_GeV\tlimit_pb\tlimit_cm2" << endl;
  const double minMass = 1.0;
  const double maxMass = 1.0e4;
  const int nPoints = 800;
  for (int i = 0; i < nPoints; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(nPoints - 1);
    const double mass = exp(log(minMass) + t * (log(maxMass) - log(minMass)));
    const double limitPb = DirDetexcl(mass);
    data << mass << "\t" << limitPb << "\t" << limitPb * 1.0e-36 << endl;
  }
  data.close();

  ofstream script(scriptPath.c_str());
  script << "set terminal pngcairo size 1000,750 enhanced font 'Arial,12'\n";
  script << "set output " << shell_quote(plotPath) << "\n";
  script << "set logscale xy\n";
  script << "set grid\n";
  script << "set xlabel 'dark matter mass m_{DM} [GeV]'\n";
  script << "set ylabel 'direct detection limit [cm^2]'\n";
  script << "set title 'DirDetexcl(m_{DM}) limit used by mO_excluder.cpp'\n";
  script << "plot " << shell_quote(dataPath)
         << " using 1:3 with lines linewidth 2 title 'DirDetexcl'\n";
  script.close();

  const string command = "gnuplot " + shell_quote(scriptPath);
  if (system(command.c_str()) != 0) {
    cerr << "Wrote " << dataPath << " and " << scriptPath
         << "; install gnuplot or run the script manually to create "
         << plotPath << endl;
  } else {
    cout << "Wrote " << dataPath << " and " << plotPath << endl;
  }
}

struct FermiLineLimitPoint {
  double energyGeV;
  double fluxLimitCm2S;
};

vector<FermiLineLimitPoint> fermi_lat_r16_line_limits() {
  const double unit = 1.0e-9;
  const double values[][2] = {
    {5.00, 3.97}, {5.20, 2.96}, {5.41, 2.25}, {5.62, 1.83},
    {5.85, 1.90}, {6.08, 2.03}, {6.33, 2.22}, {6.58, 2.21},
    {6.84, 2.06}, {7.12, 2.05}, {7.40, 1.16}, {7.70, 0.56},
    {8.01, 0.86}, {8.33, 1.51}, {8.67, 1.36}, {9.02, 1.08},
    {9.39, 0.89}, {9.77, 0.73}, {10.17, 0.51}, {10.59, 0.87},
    {11.02, 1.82}, {11.48, 1.48}, {11.96, 0.66}, {12.46, 0.50},
    {12.98, 1.26}, {13.53, 1.55}, {14.10, 0.97}, {14.70, 0.43},
    {15.33, 0.40}, {15.99, 0.37}, {16.69, 0.39}, {17.42, 0.61},
    {18.18, 0.70}, {18.99, 0.55}, {19.84, 0.46}, {20.73, 0.45},
    {21.66, 0.42}, {22.64, 0.50}, {23.66, 0.90}, {24.74, 0.77},
    {25.86, 0.81}, {27.04, 0.72}, {28.28, 0.37}, {29.57, 0.24},
    {30.93, 0.18}, {32.36, 0.19}, {33.85, 0.22}, {35.42, 0.58},
    {37.07, 0.66}, {38.80, 0.45}, {40.62, 0.36}, {42.54, 0.65},
    {44.55, 0.46}, {46.66, 0.50}, {48.88, 0.35}, {51.22, 0.16},
    {53.69, 0.21}, {56.30, 0.50}, {59.05, 0.32}, {61.96, 0.35},
    {65.04, 0.22}, {68.29, 0.37}, {71.75, 0.44}, {75.41, 0.29},
    {79.30, 0.10}, {83.43, 0.09}, {87.82, 0.18}, {92.51, 0.10},
    {97.50, 0.08}, {102.82, 0.14}, {108.49, 0.22}, {114.51, 0.33},
    {120.89, 0.42}, {127.66, 0.37}, {134.86, 0.38}, {142.51, 0.28},
    {150.66, 0.14}, {159.32, 0.18}, {168.56, 0.20}, {178.41, 0.20},
    {188.92, 0.14}, {200.15, 0.12}, {212.16, 0.14}, {225.08, 0.06},
    {239.01, 0.05}, {254.05, 0.08}, {270.33, 0.09}, {300.00, 0.13}
  };

  vector<FermiLineLimitPoint> limits;
  const int nValues = sizeof(values) / sizeof(values[0]);
  for (int i = 0; i < nValues; ++i) {
    FermiLineLimitPoint point = {values[i][0], values[i][1] * unit};
    limits.push_back(point);
  }
  return limits;
}

double fermi_lat_r16_line_limit(double energyGeV) {
  const vector<FermiLineLimitPoint> limits = fermi_lat_r16_line_limits();
  if (limits.empty() || energyGeV < limits.front().energyGeV ||
      energyGeV > limits.back().energyGeV) {
    return numeric_limits<double>::infinity();
  }

  for (size_t i = 0; i + 1 < limits.size(); ++i) {
    const FermiLineLimitPoint& low = limits[i];
    const FermiLineLimitPoint& high = limits[i + 1];
    if (energyGeV == low.energyGeV) {
      return low.fluxLimitCm2S;
    }
    if (energyGeV <= high.energyGeV) {
      const double logEnergy =
          (log(energyGeV) - log(low.energyGeV)) /
          (log(high.energyGeV) - log(low.energyGeV));
      return exp(log(low.fluxLimitCm2S) +
                 logEnergy * (log(high.fluxLimitCm2S) -
                              log(low.fluxLimitCm2S)));
    }
  }

  return limits.back().fluxLimitCm2S;
}

struct IndirectLimitResult {
  bool available;
  bool excluded;
  int channelsSeen;
  int channelsUsed;
  double maxRatio;
  double energyGeV;
  double fluxCm2S;
  double limitCm2S;
};

struct IndirectLineChannel {
  double energyGeV;
  double fluxCm2S;
};

IndirectLimitResult assess_indirect_limit(const vector<IndirectLineChannel>& channels) {
  IndirectLimitResult result;
  result.available = false;
  result.excluded = false;
  result.channelsSeen = 0;
  result.channelsUsed = 0;
  result.maxRatio = 0.0;
  result.energyGeV = numeric_limits<double>::quiet_NaN();
  result.fluxCm2S = numeric_limits<double>::quiet_NaN();
  result.limitCm2S = numeric_limits<double>::quiet_NaN();

  for (vector<IndirectLineChannel>::const_iterator it = channels.begin();
       it != channels.end(); ++it) {
    ++result.channelsSeen;

    const double energy = it->energyGeV;
    const double flux = it->fluxCm2S;
    const double limit = fermi_lat_r16_line_limit(energy);

    if (!isfinite(limit) || limit <= 0.0 || !isfinite(flux) || flux < 0.0) {
      continue;
    }

    ++result.channelsUsed;
    result.available = true;
    const double ratio = flux / limit;
    if (ratio > result.maxRatio) {
      result.maxRatio = ratio;
      result.energyGeV = energy;
      result.fluxCm2S = flux;
      result.limitCm2S = limit;
    }
  }

  result.excluded = result.available && result.maxRatio > 1.0;
  return result;
}

void write_indirect_limit_plot(const string& outputDir) {
  if (!ensure_directory(outputDir)) {
    cerr << "Could not create output directory " << outputDir << endl;
    return;
  }

  const string dataPath = output_path(outputDir, "a_fermi_lat_line_limit_r16.dat");
  const string scriptPath = output_path(outputDir, "a_fermi_lat_line_limit_r16.gnuplot");
  const string plotPath = output_path(outputDir, "a_fermi_lat_line_limit_r16.png");

  const vector<FermiLineLimitPoint> limits = fermi_lat_r16_line_limits();
  ofstream data(dataPath.c_str());
  data << "# Fermi-LAT 95% CL observed gamma-gamma line flux upper limits, R16 ROI\n";
  data << "# Source: arXiv:1305.5597, Tables VII-X. Flux column Phi_gammagamma.\n";
  data << "# E_gamma_GeV\tPhi_limit_cm-2_s-1" << endl;
  for (vector<FermiLineLimitPoint>::const_iterator it = limits.begin();
       it != limits.end(); ++it) {
    data << it->energyGeV << "\t" << it->fluxLimitCm2S << endl;
  }
  data.close();

  ofstream script(scriptPath.c_str());
  script << "set terminal pngcairo size 1000,750 enhanced font 'Arial,12'\n";
  script << "set output " << shell_quote(plotPath) << "\n";
  script << "set logscale xy\n";
  script << "set grid\n";
  script << "set xlabel 'line photon energy E_{gamma} [GeV]'\n";
  script << "set ylabel '95% CL line flux upper limit [cm^{-2} s^{-1}]'\n";
  script << "set title 'Fermi-LAT gamma-ray line limit, R16 ROI (arXiv:1305.5597)'\n";
  script << "plot " << shell_quote(dataPath)
         << " using 1:2 with linespoints linewidth 2 pointtype 7 pointsize 0.45"
         << " title 'Observed Phi_{gamma gamma} limit'\n";
  script.close();

  const string command = "gnuplot " + shell_quote(scriptPath);
  if (system(command.c_str()) != 0) {
    cerr << "Wrote " << dataPath << " and " << scriptPath
         << "; install gnuplot or run the script manually to create "
         << plotPath << endl;
  } else {
    cout << "Wrote " << dataPath << " and " << plotPath << endl;
  }
}

int main(int argc, char* argv[]){

  bool rescale=true; // multicomp dm rescaling ??
  const char* outputDirEnv = getenv("MO_EXCLUDER_OUTPUT_DIR");
  const string outputDir = outputDirEnv ? outputDirEnv : "../output";

  if (argc == 2 && string(argv[1]) == "--plot-dirdet-limits") {
    write_dirdet_limit_plot(outputDir);
    return 0;
  }

  if (argc == 2 && string(argv[1]) == "--plot-indirect-limits") {
    write_indirect_limit_plot(outputDir);
    return 0;
  }

  if (argc < 5 || ((argc - 5) % 2) != 0) {
    cout<<"Sth wrong with inp!!"<<endl;
    cout<<"i have "<<argc<<" arguments"<<endl;
    cout<<"usage: "<<argv[0]<<" index m_DM Omega DirDet [E_gamma Phi_R16]..."<<endl;
    cout<<"       "<<argv[0]<<" --plot-dirdet-limits"<<endl;
    cout<<"       "<<argv[0]<<" --plot-indirect-limits"<<endl;
    return 1;
  }

  int index= (int)atof(argv[1]);
  double m_DM  = (double)atof(argv[2]);
  double Omega = (double)atof(argv[3]);
  double DirDet = (double)atof(argv[4]);

  vector<IndirectLineChannel> indirectChannels;
  for (int i = 5; i + 1 < argc; i += 2) {
    IndirectLineChannel channel = {atof(argv[i]), atof(argv[i + 1])};
    indirectChannels.push_back(channel);
  }

  cout<<"processing "<<index<<"\t"<<m_DM<<"\t"<<Omega<<"\t"<<DirDet<<endl;

  double lx = 0.0, lhx = 0.0, lsx = 0.0, mx = 0.0, vevs = 0.0, sint = 0.0, mh2 = 0.0;
  bool foundPoint = false;
  const char* oksFileEnv = getenv("MO_EXCLUDER_OKS_FILE");
  const string oksFile = oksFileEnv ? oksFileEnv : "../run/oks.dat";

  ifstream oksfile(oksFile.c_str());
  if (!oksfile) {
    cerr << "Could not open " << oksFile << endl;
    return 1;
  }

  int ii = 0;
  double lxIn, lhxIn, lsxIn, mxIn, vevsIn, sintIn, mh2In;
  while (oksfile >> ii >> lxIn >> lhxIn >> lsxIn >> mxIn >> vevsIn >> sintIn >> mh2In) {
    if (ii == index) {
      lx = lxIn;
      lhx = lhxIn;
      lsx = lsxIn;
      mx = mxIn;
      vevs = vevsIn;
      sint = sintIn;
      mh2 = mh2In;
      foundPoint = true;
      break;
    }
  }
  oksfile.close();

  if (!foundPoint) {
    cerr << "Index " << index << " not found in oks.dat" << endl;
    return 1;
  }

  if (!ensure_directory(outputDir)) {
    cerr << "Could not create output directory " << outputDir << endl;
    return 1;
  }

  bool DMisok = true;
  bool OMGok = true;
  bool DIRok = true;
  bool INDok = true;

  const double relicUpperLimit = 0.121;
  const double relicStrictCentral = 0.12;
  const double relicStrictWidth = 0.001;
  const double LUX = DirDetexcl(m_DM);
  double dirDetLimit = LUX;

  if (rescale == true) {
    if (Omega > 0.0) {
      dirDetLimit = LUX * relicUpperLimit / Omega;
    } else {
      dirDetLimit = numeric_limits<double>::infinity();
    }
  }

  ostringstream baseLine;
  baseLine << index << "\t" << lx << "\t" << lhx << "\t" << lsx << "\t" << mx
           << "\t" << vevs << "\t" << sint << "\t" << mh2 << "\t" << m_DM;

  const IndirectLimitResult indirect = assess_indirect_limit(indirectChannels);

  ostringstream dmLine;
  dmLine << baseLine.str() << "\t" << Omega << "\t" << DirDet << "\t"
         << dirDetLimit << "\t" << LUX;

  ostringstream indirectLine;
  indirectLine << dmLine.str() << "\t" << (indirect.available ? 1 : 0)
               << "\t" << indirect.energyGeV << "\t" << indirect.fluxCm2S
               << "\t" << indirect.limitCm2S << "\t" << indirect.maxRatio;

  ofstream DM(output_path(outputDir + "/DM_data", "DM_data").c_str());
  DM << dmLine.str() << endl;
  DM.close();

  if (Omega <= relicUpperLimit) {
    append_line(output_path(outputDir, "relic_pass.dat"), dmLine.str());
  }

  if (fabs(Omega - relicStrictCentral) <= relicStrictWidth) {
    append_line(output_path(outputDir, "relic_strict.dat"), dmLine.str());
  }

  if (Omega > relicUpperLimit) {
    DMisok = false;
    OMGok = false;
    ostringstream omLine;
    omLine << baseLine.str() << "\t" << Omega;
    append_line(output_path(outputDir, "omexcl.dat"), omLine.str());
  }

  if (DirDet > dirDetLimit) {
    DMisok = false;
    DIRok = false;
    ostringstream luxLine;
    luxLine << baseLine.str() << "\t" << Omega << "\t" << DirDet << "\t"
            << dirDetLimit << "\t" << LUX;
    append_line(output_path(outputDir, "luxexcl.dat"), luxLine.str());
  }

  if (OMGok && DIRok) {
    append_line(output_path(outputDir, "all_dirpass.dat"), dmLine.str());
  }

  if (indirect.excluded) {
    DMisok = false;
    INDok = false;
    append_line(output_path(outputDir, "indirexcl.dat"), indirectLine.str());
  } else {
    append_line(output_path(outputDir, "indirpass.dat"), indirectLine.str());
  }

  if (!DMisok) {
    cout<<"excluded from dm "<<OMGok<<DIRok<<INDok<<endl;
    append_line(output_path(outputDir, "dmexcl.dat"), dmLine.str());

    // ofstream excl(output_path(outputDir + "/DM_EXCLUDED", "DM_EXCLUDED").c_str());
    // excl.close();

    if (!OMGok) {
      // ofstream RelDens(output_path(outputDir + "/RelDens_EXCLUDED", "RelDens_EXCLUDED").c_str());
      // RelDens.close();
      cout<<"Too high relic dens"<<endl;
    }

    if (!DIRok) {
      // ofstream DirDetExcl(output_path(outputDir + "/DirDet_EXCLUDED", "DirDet_EXCLUDED").c_str());
      // DirDetExcl.close();
      cout<<"Should be visible in LUX!"<<endl;
    }

    if (!INDok) {
      // ofstream IndirDetExcl(output_path(outputDir + "/IndirDet_EXCLUDED", "IndirDet_EXCLUDED").c_str());
      // IndirDetExcl.close();
      cout<<"Excluded by Fermi-LAT gamma-line indirect detection"<<endl;
    }
  } else {
    cout<<"Point agrees with DM data"<<endl;
    append_line(output_path(outputDir, "allall.dat"), dmLine.str());
  }

  return 0;

}
