#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <cerrno>
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

}  // namespace


/*
double LUXexcl(double mH){
double excl;
// for the fitted function see LUX_fit
excl = 1.86934E-10 + 1.1236E-11*mH + (-1.32312E-16)*mH*mH + (1.85167E8)/((-9.35439E15)*mH + 1.12151E15*mH*mH);
return excl;} 

*/

double LUXexcl(double mH){                                                                                                                                                 double excl; 

  if ((mH < 17.) || mH > 200.)
    excl=-0.138541 + 0.00343701*mH - 1.26198e-6*mH*mH + 3.44392e+10/(-1.49313e+9*mH + 2.32828e+8*mH*mH);
  else if (mH < 33.)
    excl=-0.0209648 - 0.00479147*mH + 0.0000388817*mH*mH + 2.75483/(0.248816*mH - 0.00173454*mH*mH);
  else excl=0.295076 - 0.00150593*mH + 0.0000109119*mH*mH + 0.999075/(1.00003*mH + 1.0009*mH*mH);
  excl=excl*1.e-9;
  return excl;}

double XENONexcl(double x){

  double result;
  if (x<10){
    result= (exp(-23.2949 + 2.2493* pow((-3.60228 + log(x)),2)));
    result=result*(1-0.465188 + 1.64349*pow((-2.3605 + log(x)),2));}
  // else if (x<11) result=(8.75459*1.e-9 - 6.96272*1.e-10* x);
  //  else if (x<18) {result=exp(-23.2938 + 0.0135204* pow((-24.8742 + x),2));
  //  result=result*(1-0.0826054 + 0.0110134*pow((-13.3867 + x),2));}
  // else if (x<20)result= 5.0606*1.e-10 - 1.91361*1.e-11* x;
  else if (x<20) result=(9.301626524677054e-8*exp(33.018821405901335/log(x))*pow(x,1.5351327504595669))*1.e-10;
  // else if (x<30) result=5.103100000000005e-10 - 3.674166666666671e-11*x + 7.027333333333344e-13*pow(x,2);
  else if (x<30) result=1.744393487921861e-28*exp(68.14955833703276/log(x))*pow(x,5.8655281320921215);

  //  else if (x<30) result=(7.64124*1.e-11 + 1.34496*1.e-13*pow ((-36.6817 + x),2));


  else if (x<40) result= (2.771+0.0427*x)*1.e-11;
  //  else if (x<65) result=3.050728*1.e-13*x-3.13092*1.e-12; // new lux zeplin for this region only
    else if (x< 60) result= 4.84775*1.e-11 + 7.28365*1.e-13* x;

  else if (x<100) result = exp(0.671326*log(x)-25.8196);
  else result=1.49167*1.e-12*exp(0.977674*log(x));

  //  if (( x> 400) && (result > 8.243e-10)) result=8.243e-10; //very rough cutoff for new xenon1t

  if ((x> 40) && (x<100)) result=(0.0774*x+1.382)*1.e-11;
  
  if ((x> 50) && (x<70)) result=(0.0767*x+1.423)*1.e-11;
  else if ((x > 100) && (x <200)) result = (0.0805*x+1.235)*1.e-11; 
  else if ((x>200) && (x<400)) result=8.135e-13*x+8.8e-12;
  else if (x>400) result=8.168e-13*x+7.5e-12; // very rough approx for new xenon1t

  if ((x>40) && (x<65)) result=3.050728*1.e-13*x-3.13092*1.e-12; // new lux zeplin for this region only    
  if ((x>= 65) && (x<110)) result=2.80493*1.e-13*x-1.53327*1.e-12; // same
  if ((x>=78) && (x < 219)) result=(0.0281889*x-0.154847)*1.e-11; // same
  
  return result;


  //  double result = exp(0.671326*log(mH)-25.8196);
  // return result;
}



/*


double XENONexcl(double x){

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


int main(int argc, char* argv[]){

  bool rescale=true; // multicomp dm rescaling ??
  if(argc!=5){
    cout<<"Sth wrong with inp!!"<<endl;
    cout<<"i have "<<argc<<" arguments"<<endl;
    return 1;
  }

  int index= (int)atof(argv[1]);
  double m_DM  = (double)atof(argv[2]);
  double Omega = (double)atof(argv[3]);
  double DirDet = (double)atof(argv[4]);

  cout<<"processing "<<index<<"\t"<<m_DM<<"\t"<<Omega<<"\t"<<DirDet<<endl;

  double lx = 0.0, lhx = 0.0, lsx = 0.0, mx = 0.0, vevs = 0.0, sint = 0.0, mh2 = 0.0;
  bool foundPoint = false;
  const char* outputDirEnv = getenv("MO_EXCLUDER_OUTPUT_DIR");
  const string outputDir = outputDirEnv ? outputDirEnv : "../output";
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

  const double relicUpperLimit = 0.1224;
  const double LUX = XENONexcl(m_DM);
  double dirDetLimit = LUX;

  if (rescale == true) {
    if (Omega > 0.0) {
      dirDetLimit = LUX * relicUpperLimit / Omega;
    } else {
      dirDetLimit = numeric_limits<double>::infinity();
    }
  }

  ofstream DM(output_path(outputDir, "DM_data").c_str());
  DM << index << "\t" << lx << "\t" << lhx << "\t" << lsx << "\t" << mx << "\t"
     << vevs << "\t" << sint << "\t" << mh2 << "\t" << m_DM << "\t" << Omega
     << "\t" << DirDet << "\t" << dirDetLimit << "\t" << LUX << endl;
  DM.close();

  ostringstream baseLine;
  baseLine << index << "\t" << lx << "\t" << lhx << "\t" << lsx << "\t" << mx
           << "\t" << vevs << "\t" << sint << "\t" << mh2 << "\t" << m_DM;

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

  if (!DMisok) {
    cout<<"excluded from dm "<<OMGok<<DIRok<<endl;
    ostringstream dmLine;
    dmLine << baseLine.str() << "\t" << Omega << "\t" << DirDet << "\t"
           << dirDetLimit << "\t" << LUX;
    append_line(output_path(outputDir, "dmexcl.dat"), dmLine.str());

    ofstream excl(output_path(outputDir, "DM_EXCLUDED").c_str());
    excl.close();

    if (!OMGok) {
      ofstream RelDens(output_path(outputDir, "RelDens_EXCLUDED").c_str());
      RelDens.close();
      cout<<"Too high relic dens"<<endl;
    }

    if (!DIRok) {
      ofstream DirDetExcl(output_path(outputDir, "DirDet_EXCLUDED").c_str());
      DirDetExcl.close();
      cout<<"Should be visible in LUX!"<<endl;
    }
  } else {
    cout<<"Point agrees with DM data"<<endl;
    ostringstream okLine;
    okLine << baseLine.str() << "\t" << Omega << "\t" << DirDet << "\t"
           << dirDetLimit << "\t" << LUX;
    append_line(output_path(outputDir, "allall.dat"), okLine.str());
  }

  return 0;

}
