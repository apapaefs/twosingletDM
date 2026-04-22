#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

namespace {

bool ensure_directory(const string& path) {
  if (mkdir(path.c_str(), 0777) == 0) {
    return true;
  }

  return errno == EEXIST;
}

}  // namespace

int main(void) {
  ifstream pointsin("../run/oks.dat", ios::in);
  if (!pointsin) {
    cerr << "Failed to open ../run/oks.dat" << endl;
    return 1;
  }

  if (!ensure_directory("../run") || !ensure_directory("../run/cards")) {
    cerr << "Failed to create ../run/cards" << endl;
    return 1;
  }

  int ii;
  double lx_in, lhx_in, lsx_in, mx_in, vevs_in, sint_in, mh2_in;

  while (pointsin >> ii >> lx_in >> lhx_in >> lsx_in >> mx_in >> vevs_in >>
         sint_in >> mh2_in) {
    cout << ii << endl;

    stringstream convert;
    convert << ii;
    const string filepattern = convert.str() + ".dat";

    ofstream mO_inp("../run/cards/MO_inp" + filepattern);
    if (!mO_inp) {
      cerr << "Failed to open ../run/cards/MO_inp" << filepattern << endl;
      return 1;
    }

    mO_inp << "LX\t\t" << lx_in << endl;
    mO_inp << "LHX\t\t" << lhx_in << endl;
    mO_inp << "LSX\t\t" << lsx_in << endl;
    mO_inp << "MX\t\t" << mx_in << endl;
    mO_inp << "vevs\t" << vevs_in << endl;
    mO_inp << "SinT\t" << sint_in << endl;
    mO_inp << "Mh2\t\t" << mh2_in << endl;
  }

  return 0;
}
