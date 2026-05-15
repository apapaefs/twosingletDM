#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialTRSM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_TRSM::Class_Potential_TRSM(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::TRSM;

  nPar = 9;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 15; // number of parameters in the counterterm potential

  nVEV = 3; // number of VEVs to minimize the potential

  NHiggs = 6; // number of scalar d.o.f.

  NGauge = 4; // number of gauge fields

  NLepton = 9; // number of lepton fields

  NQuarks = 12; // number of quark fields

  VevOrder.resize(nVEV);
  VevOrder[0] = 3; // w1
  VevOrder[1] = 4; // wx
  VevOrder[2] = 5; // ws

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_TRSM::~Class_Potential_TRSM()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_TRSM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmphiSq");
  labels.push_back("dlphi");
  labels.push_back("dmsSq");
  labels.push_back("dls");
  labels.push_back("dmxSq");
  labels.push_back("dlx");
  labels.push_back("dlphis");
  labels.push_back("dlphix");
  labels.push_back("dlsx");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");
  labels.push_back("dT6");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_TRSM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");     // Label for the critical temperature
  labels.push_back("v_c");     // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for xi_c
  // out += "VEV order";
  labels.push_back("w1(T_c)");
  labels.push_back("wx(T_c)");
  labels.push_back("ws(T_c)");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_TRSM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles
  particles.push_back("h_1");
  particles.push_back("h_2");
  particles.push_back("h_3");
  particles.push_back("h_4");
  particles.push_back("h_5");
  particles.push_back("h_6");

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_TRSM::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order"
  labels.push_back("w1");
  labels.push_back("wx");
  labels.push_back("ws");

  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_TRSM::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 8; k++)
  {
    ss >> tmp;
    if (k == 1)
      par[0] = tmp; // m1
    if (k == 2)
      par[1] = tmp; // m2
    if (k == 3)
      par[2] = tmp; // m3
    if (k == 4)
      par[3] = tmp; // vs
    if (k == 5)
      par[4] = tmp; // a12
    if (k == 6)
      par[5] = tmp; // lx
    if (k == 7)
      par[6] = tmp; // lphix
    if (k == 8)
      par[7] = tmp; // lsx

  }

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_TRSM::set_gen(const std::vector<double> &par)
{

  m1 = par[0]; 
  m2 = par[1]; 
  m3 = par[2]; 
  vs = par[3]; 
  a12 = par[4]; 
  lx = par[5]; 
  lphix = par[6]; 
  lsx = par[7]; 

  v1 = SMConstants.C_vev0; 
  lphi = (pow(m1,2)*pow(cos(a12),2) + pow(m2,2)*pow(sin(a12),2))/(2.*pow(SMConstants.C_vev0,2)); 
  ls = (pow(m2,2)*pow(cos(a12),2) + pow(m1,2)*pow(sin(a12),2))/(8.*pow(vs,2)); 
  lphis = ((-pow(m1,2) + pow(m2,2))*cos(a12)*sin(a12))/(2.*SMConstants.C_vev0*vs); 
  mxSq = (pow(m3,2) - lphix*pow(SMConstants.C_vev0,2) - 2*lsx*pow(vs,2))/2.; 
  mphiSq = -0.5*(pow(m1,2)*SMConstants.C_vev0*pow(cos(a12),2) - (pow(m1,2) - pow(m2,2))*vs*cos(a12)*sin(a12) + pow(m2,2)*SMConstants.C_vev0*pow(sin(a12),2))/SMConstants.C_vev0; 
  msSq = -0.25*(pow(m2,2)*vs*pow(cos(a12),2) - (pow(m1,2) - pow(m2,2))*SMConstants.C_vev0*cos(a12)*sin(a12) + pow(m1,2)*vs*pow(sin(a12),2))/vs; 

  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV
  vevTreeMin[0] = v1; // w1
  vevTreeMin[1] = 0; // wx
  vevTreeMin[2] = vs; // ws

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_TRSM::set_CT_Pot_Par(const std::vector<double> &par)
{
  dmphiSq = par[0];
  dlphi = par[1];
  dmsSq = par[2];
  dls = par[3];
  dmxSq = par[4];
  dlx = par[5];
  dlphis = par[6];
  dlphix = par[7];
  dlsx = par[8];
  dT1 = par[9];
  dT2 = par[10];
  dT3 = par[11];
  dT4 = par[12];
  dT5 = par[13];
  dT6 = par[14];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L1[5] = dT6;

  Curvature_Higgs_CT_L2[0][0] = dmphiSq;
  Curvature_Higgs_CT_L2[1][1] = dmphiSq;
  Curvature_Higgs_CT_L2[2][2] = dmphiSq;
  Curvature_Higgs_CT_L2[3][3] = dmphiSq;
  Curvature_Higgs_CT_L2[4][4] = 2*dmxSq;
  Curvature_Higgs_CT_L2[5][5] = 2*dmsSq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 6*dlphi;
  Curvature_Higgs_CT_L4[0][0][1][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][0][2][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][0][3][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][0][4][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[0][0][5][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[0][1][0][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][1][1][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][2][0][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][2][2][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][3][0][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][3][3][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[0][4][0][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[0][4][4][0] = 2*dlphix;
  Curvature_Higgs_CT_L4[0][5][0][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[0][5][5][0] = 2*dlphis;
  Curvature_Higgs_CT_L4[1][0][0][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][0][1][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][1][0][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][1][1][1] = 6*dlphi;
  Curvature_Higgs_CT_L4[1][1][2][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][1][3][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][1][4][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[1][1][5][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[1][2][1][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][2][2][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][3][1][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][3][3][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[1][4][1][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[1][4][4][1] = 2*dlphix;
  Curvature_Higgs_CT_L4[1][5][1][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[1][5][5][1] = 2*dlphis;
  Curvature_Higgs_CT_L4[2][0][0][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][0][2][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][1][1][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][1][2][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][2][0][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][2][1][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][2][2][2] = 6*dlphi;
  Curvature_Higgs_CT_L4[2][2][3][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][2][4][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[2][2][5][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[2][3][2][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][3][3][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[2][4][2][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[2][4][4][2] = 2*dlphix;
  Curvature_Higgs_CT_L4[2][5][2][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[2][5][5][2] = 2*dlphis;
  Curvature_Higgs_CT_L4[3][0][0][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][0][3][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][1][1][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][1][3][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][2][2][3] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][2][3][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][3][0][0] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][3][1][1] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][3][2][2] = 2*dlphi;
  Curvature_Higgs_CT_L4[3][3][3][3] = 6*dlphi;
  Curvature_Higgs_CT_L4[3][3][4][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[3][3][5][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[3][4][3][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[3][4][4][3] = 2*dlphix;
  Curvature_Higgs_CT_L4[3][5][3][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[3][5][5][3] = 2*dlphis;
  Curvature_Higgs_CT_L4[4][0][0][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][0][4][0] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][1][1][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][1][4][1] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][2][2][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][2][4][2] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][3][3][4] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][3][4][3] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][4][0][0] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][4][1][1] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][4][2][2] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][4][3][3] = 2*dlphix;
  Curvature_Higgs_CT_L4[4][4][4][4] = 24*dlx;
  Curvature_Higgs_CT_L4[4][4][5][5] = 4*dlsx;
  Curvature_Higgs_CT_L4[4][5][4][5] = 4*dlsx;
  Curvature_Higgs_CT_L4[4][5][5][4] = 4*dlsx;
  Curvature_Higgs_CT_L4[5][0][0][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][0][5][0] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][1][1][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][1][5][1] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][2][2][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][2][5][2] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][3][3][5] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][3][5][3] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][4][4][5] = 4*dlsx;
  Curvature_Higgs_CT_L4[5][4][5][4] = 4*dlsx;
  Curvature_Higgs_CT_L4[5][5][0][0] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][5][1][1] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][5][2][2] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][5][3][3] = 2*dlphis;
  Curvature_Higgs_CT_L4[5][5][4][4] = 4*dlsx;
  Curvature_Higgs_CT_L4[5][5][5][5] = 24*dls;


}

void Class_Potential_TRSM::AdjustRotationMatrix()
{
}

/**
 * console output of all parameters
 */
void Class_Potential_TRSM::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << "\n";

  ss << "\nThe input parameters are : \n";
  ss << "m1 = " << m1 << "\n";
  ss << "m2 = " << m2 << "\n";
  ss << "m3 = " << m3 << "\n";
  ss << "vs = " << vs << "\n";
  ss << "a12 = " << a12 << "\n";
  ss << "lx = " << lx << "\n";
  ss << "lphix = " << lphix << "\n";
  ss << "lsx = " << lsx << "\n";

  ss << "\nThe parameters are : \n";
  ss << "mphiSq = " << mphiSq << "\n";
  ss << "lphi = " << lphi << "\n";
  ss << "msSq = " << msSq << "\n";
  ss << "ls = " << ls << "\n";
  ss << "mxSq = " << mxSq << "\n";
  ss << "lx = " << lx << "\n";
  ss << "lphis = " << lphis << "\n";
  ss << "lphix = " << lphix << "\n";
  ss << "lsx = " << lsx << "\n";

  ss << "\nThe counterterm parameters are : \n";
  ss << "dmphiSq = " << dmphiSq << "\n";
  ss << "dlphi = " << dlphi << "\n";
  ss << "dmsSq = " << dmsSq << "\n";
  ss << "dls = " << dls << "\n";
  ss << "dmxSq = " << dmxSq << "\n";
  ss << "dlx = " << dlx << "\n";
  ss << "dlphis = " << dlphis << "\n";
  ss << "dlphix = " << dlphix << "\n";
  ss << "dlsx = " << dlsx << "\n";
  ss << "dT1 = " << dT1 << "\n";
  ss << "dT2 = " << dT2 << "\n";
  ss << "dT3 = " << dT3 << "\n";
  ss << "dT4 = " << dT4 << "\n";
  ss << "dT5 = " << dT5 << "\n";
  ss << "dT6 = " << dT6 << "\n";

  ss << "\nThe scale is given by mu = " << scale << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_TRSM::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // formulae for the counterterm scheme
  parCT.push_back((v1*(-3*HesseWeinberg(0,0) + HesseWeinberg(3,3)) + vs*HesseWeinberg(3,5))/(2.*v1)); //dmphiSq;
  parCT.push_back((HesseWeinberg(0,0) - HesseWeinberg(3,3))/(2.*pow(v1,2))); //dlphi;
  parCT.push_back((v1*HesseWeinberg(3,5) + vs*HesseWeinberg(5,5) - 3*NablaWeinberg(5))/(4.*vs)); //dmsSq;
  parCT.push_back((-(vs*HesseWeinberg(5,5)) + NablaWeinberg(5))/(8.*pow(vs,3))); //dls;
  parCT.push_back(-0.5*HesseWeinberg(4,4)); //dmxSq;
  parCT.push_back(0); //dlx;
  parCT.push_back(-0.5*HesseWeinberg(3,5)/(v1*vs)); //dlphis;
  parCT.push_back(0); //dlphix;
  parCT.push_back(0); //dlsx;
  parCT.push_back(-NablaWeinberg(0)); //dT1;
  parCT.push_back(-NablaWeinberg(1)); //dT2;
  parCT.push_back(-NablaWeinberg(2)); //dT3;
  parCT.push_back(v1*HesseWeinberg(0,0) - NablaWeinberg(3)); //dT4;
  parCT.push_back(-NablaWeinberg(4)); //dT5;
  parCT.push_back(0); //dT6;

  return parCT;
}

// mass basis triple couplings
void Class_Potential_TRSM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  // new rotation matrix with
  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  std::vector<double> HiggsOrder(NHiggs);

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_Potential_TRSM::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = mphiSq;
  Curvature_Higgs_L2[1][1] = mphiSq;
  Curvature_Higgs_L2[2][2] = mphiSq;
  Curvature_Higgs_L2[3][3] = mphiSq;
  Curvature_Higgs_L2[4][4] = 2*mxSq;
  Curvature_Higgs_L2[5][5] = 2*msSq;

  Curvature_Higgs_L4[0][0][0][0] = 6*lphi;
  Curvature_Higgs_L4[0][0][1][1] = 2*lphi;
  Curvature_Higgs_L4[0][0][2][2] = 2*lphi;
  Curvature_Higgs_L4[0][0][3][3] = 2*lphi;
  Curvature_Higgs_L4[0][0][4][4] = 2*lphix;
  Curvature_Higgs_L4[0][0][5][5] = 2*lphis;
  Curvature_Higgs_L4[0][1][0][1] = 2*lphi;
  Curvature_Higgs_L4[0][1][1][0] = 2*lphi;
  Curvature_Higgs_L4[0][2][0][2] = 2*lphi;
  Curvature_Higgs_L4[0][2][2][0] = 2*lphi;
  Curvature_Higgs_L4[0][3][0][3] = 2*lphi;
  Curvature_Higgs_L4[0][3][3][0] = 2*lphi;
  Curvature_Higgs_L4[0][4][0][4] = 2*lphix;
  Curvature_Higgs_L4[0][4][4][0] = 2*lphix;
  Curvature_Higgs_L4[0][5][0][5] = 2*lphis;
  Curvature_Higgs_L4[0][5][5][0] = 2*lphis;
  Curvature_Higgs_L4[1][0][0][1] = 2*lphi;
  Curvature_Higgs_L4[1][0][1][0] = 2*lphi;
  Curvature_Higgs_L4[1][1][0][0] = 2*lphi;
  Curvature_Higgs_L4[1][1][1][1] = 6*lphi;
  Curvature_Higgs_L4[1][1][2][2] = 2*lphi;
  Curvature_Higgs_L4[1][1][3][3] = 2*lphi;
  Curvature_Higgs_L4[1][1][4][4] = 2*lphix;
  Curvature_Higgs_L4[1][1][5][5] = 2*lphis;
  Curvature_Higgs_L4[1][2][1][2] = 2*lphi;
  Curvature_Higgs_L4[1][2][2][1] = 2*lphi;
  Curvature_Higgs_L4[1][3][1][3] = 2*lphi;
  Curvature_Higgs_L4[1][3][3][1] = 2*lphi;
  Curvature_Higgs_L4[1][4][1][4] = 2*lphix;
  Curvature_Higgs_L4[1][4][4][1] = 2*lphix;
  Curvature_Higgs_L4[1][5][1][5] = 2*lphis;
  Curvature_Higgs_L4[1][5][5][1] = 2*lphis;
  Curvature_Higgs_L4[2][0][0][2] = 2*lphi;
  Curvature_Higgs_L4[2][0][2][0] = 2*lphi;
  Curvature_Higgs_L4[2][1][1][2] = 2*lphi;
  Curvature_Higgs_L4[2][1][2][1] = 2*lphi;
  Curvature_Higgs_L4[2][2][0][0] = 2*lphi;
  Curvature_Higgs_L4[2][2][1][1] = 2*lphi;
  Curvature_Higgs_L4[2][2][2][2] = 6*lphi;
  Curvature_Higgs_L4[2][2][3][3] = 2*lphi;
  Curvature_Higgs_L4[2][2][4][4] = 2*lphix;
  Curvature_Higgs_L4[2][2][5][5] = 2*lphis;
  Curvature_Higgs_L4[2][3][2][3] = 2*lphi;
  Curvature_Higgs_L4[2][3][3][2] = 2*lphi;
  Curvature_Higgs_L4[2][4][2][4] = 2*lphix;
  Curvature_Higgs_L4[2][4][4][2] = 2*lphix;
  Curvature_Higgs_L4[2][5][2][5] = 2*lphis;
  Curvature_Higgs_L4[2][5][5][2] = 2*lphis;
  Curvature_Higgs_L4[3][0][0][3] = 2*lphi;
  Curvature_Higgs_L4[3][0][3][0] = 2*lphi;
  Curvature_Higgs_L4[3][1][1][3] = 2*lphi;
  Curvature_Higgs_L4[3][1][3][1] = 2*lphi;
  Curvature_Higgs_L4[3][2][2][3] = 2*lphi;
  Curvature_Higgs_L4[3][2][3][2] = 2*lphi;
  Curvature_Higgs_L4[3][3][0][0] = 2*lphi;
  Curvature_Higgs_L4[3][3][1][1] = 2*lphi;
  Curvature_Higgs_L4[3][3][2][2] = 2*lphi;
  Curvature_Higgs_L4[3][3][3][3] = 6*lphi;
  Curvature_Higgs_L4[3][3][4][4] = 2*lphix;
  Curvature_Higgs_L4[3][3][5][5] = 2*lphis;
  Curvature_Higgs_L4[3][4][3][4] = 2*lphix;
  Curvature_Higgs_L4[3][4][4][3] = 2*lphix;
  Curvature_Higgs_L4[3][5][3][5] = 2*lphis;
  Curvature_Higgs_L4[3][5][5][3] = 2*lphis;
  Curvature_Higgs_L4[4][0][0][4] = 2*lphix;
  Curvature_Higgs_L4[4][0][4][0] = 2*lphix;
  Curvature_Higgs_L4[4][1][1][4] = 2*lphix;
  Curvature_Higgs_L4[4][1][4][1] = 2*lphix;
  Curvature_Higgs_L4[4][2][2][4] = 2*lphix;
  Curvature_Higgs_L4[4][2][4][2] = 2*lphix;
  Curvature_Higgs_L4[4][3][3][4] = 2*lphix;
  Curvature_Higgs_L4[4][3][4][3] = 2*lphix;
  Curvature_Higgs_L4[4][4][0][0] = 2*lphix;
  Curvature_Higgs_L4[4][4][1][1] = 2*lphix;
  Curvature_Higgs_L4[4][4][2][2] = 2*lphix;
  Curvature_Higgs_L4[4][4][3][3] = 2*lphix;
  Curvature_Higgs_L4[4][4][4][4] = 24*lx;
  Curvature_Higgs_L4[4][4][5][5] = 4*lsx;
  Curvature_Higgs_L4[4][5][4][5] = 4*lsx;
  Curvature_Higgs_L4[4][5][5][4] = 4*lsx;
  Curvature_Higgs_L4[5][0][0][5] = 2*lphis;
  Curvature_Higgs_L4[5][0][5][0] = 2*lphis;
  Curvature_Higgs_L4[5][1][1][5] = 2*lphis;
  Curvature_Higgs_L4[5][1][5][1] = 2*lphis;
  Curvature_Higgs_L4[5][2][2][5] = 2*lphis;
  Curvature_Higgs_L4[5][2][5][2] = 2*lphis;
  Curvature_Higgs_L4[5][3][3][5] = 2*lphis;
  Curvature_Higgs_L4[5][3][5][3] = 2*lphis;
  Curvature_Higgs_L4[5][4][4][5] = 4*lsx;
  Curvature_Higgs_L4[5][4][5][4] = 4*lsx;
  Curvature_Higgs_L4[5][5][0][0] = 2*lphis;
  Curvature_Higgs_L4[5][5][1][1] = 2*lphis;
  Curvature_Higgs_L4[5][5][2][2] = 2*lphis;
  Curvature_Higgs_L4[5][5][3][3] = 2*lphis;
  Curvature_Higgs_L4[5][5][4][4] = 4*lsx;
  Curvature_Higgs_L4[5][5][5][5] = 24*ls;

  Curvature_Gauge_G2H2[0][0][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][3][0][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][1][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][2][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][3][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][1][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][3][0][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][1][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][2][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][3][1] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][2][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][3][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][2][2] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][3][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][0][0][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][1][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][2][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][3][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][0][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][1][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][2][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][3][1] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][2][2] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][3][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][3][0][0] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][1][1] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][2][2] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][3][3] = pow(SMConstants.C_gs,2)/2.;

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  Curvature_Lepton_F2H1[0][1][2] = (II*SMConstants.C_MassElectron)/v1;
  Curvature_Lepton_F2H1[0][1][3] = SMConstants.C_MassElectron/v1;
  Curvature_Lepton_F2H1[1][0][2] = (II*SMConstants.C_MassElectron)/v1;
  Curvature_Lepton_F2H1[1][0][3] = SMConstants.C_MassElectron/v1;
  Curvature_Lepton_F2H1[1][6][0] = SMConstants.C_MassElectron/v1;
  Curvature_Lepton_F2H1[1][6][1] = (II*SMConstants.C_MassElectron)/v1;
  Curvature_Lepton_F2H1[2][3][2] = (II*SMConstants.C_MassMu)/v1;
  Curvature_Lepton_F2H1[2][3][3] = SMConstants.C_MassMu/v1;
  Curvature_Lepton_F2H1[3][2][2] = (II*SMConstants.C_MassMu)/v1;
  Curvature_Lepton_F2H1[3][2][3] = SMConstants.C_MassMu/v1;
  Curvature_Lepton_F2H1[3][7][0] = SMConstants.C_MassMu/v1;
  Curvature_Lepton_F2H1[3][7][1] = (II*SMConstants.C_MassMu)/v1;
  Curvature_Lepton_F2H1[4][5][2] = (II*SMConstants.C_MassTau)/v1;
  Curvature_Lepton_F2H1[4][5][3] = SMConstants.C_MassTau/v1;
  Curvature_Lepton_F2H1[5][4][2] = (II*SMConstants.C_MassTau)/v1;
  Curvature_Lepton_F2H1[5][4][3] = SMConstants.C_MassTau/v1;
  Curvature_Lepton_F2H1[5][8][0] = SMConstants.C_MassTau/v1;
  Curvature_Lepton_F2H1[5][8][1] = (II*SMConstants.C_MassTau)/v1;
  Curvature_Lepton_F2H1[6][1][0] = SMConstants.C_MassElectron/v1;
  Curvature_Lepton_F2H1[6][1][1] = (II*SMConstants.C_MassElectron)/v1;
  Curvature_Lepton_F2H1[7][3][0] = SMConstants.C_MassMu/v1;
  Curvature_Lepton_F2H1[7][3][1] = (II*SMConstants.C_MassMu)/v1;
  Curvature_Lepton_F2H1[8][5][0] = SMConstants.C_MassTau/v1;
  Curvature_Lepton_F2H1[8][5][1] = (II*SMConstants.C_MassTau)/v1;

  Curvature_Quark_F2H1[0][6][2] = (-II*SMConstants.C_MassUp)/v1;
  Curvature_Quark_F2H1[0][6][3] = SMConstants.C_MassUp/v1;
  Curvature_Quark_F2H1[0][9][0] = -((SMConstants.C_MassUp*conj(V11))/v1);
  Curvature_Quark_F2H1[0][9][1] = (II*SMConstants.C_MassUp*conj(V11))/v1;
  Curvature_Quark_F2H1[0][10][0] = -((SMConstants.C_MassUp*conj(V12))/v1);
  Curvature_Quark_F2H1[0][10][1] = (II*SMConstants.C_MassUp*conj(V12))/v1;
  Curvature_Quark_F2H1[0][11][0] = -((SMConstants.C_MassUp*conj(V13))/v1);
  Curvature_Quark_F2H1[0][11][1] = (II*SMConstants.C_MassUp*conj(V13))/v1;
  Curvature_Quark_F2H1[1][7][2] = (-II*SMConstants.C_MassCharm)/v1;
  Curvature_Quark_F2H1[1][7][3] = SMConstants.C_MassCharm/v1;
  Curvature_Quark_F2H1[1][9][0] = -((SMConstants.C_MassCharm*conj(V21))/v1);
  Curvature_Quark_F2H1[1][9][1] = (II*SMConstants.C_MassCharm*conj(V21))/v1;
  Curvature_Quark_F2H1[1][10][0] = -((SMConstants.C_MassCharm*conj(V22))/v1);
  Curvature_Quark_F2H1[1][10][1] = (II*SMConstants.C_MassCharm*conj(V22))/v1;
  Curvature_Quark_F2H1[1][11][0] = -((SMConstants.C_MassCharm*conj(V23))/v1);
  Curvature_Quark_F2H1[1][11][1] = (II*SMConstants.C_MassCharm*conj(V23))/v1;
  Curvature_Quark_F2H1[2][8][2] = (-II*SMConstants.C_MassTop)/v1;
  Curvature_Quark_F2H1[2][8][3] = SMConstants.C_MassTop/v1;
  Curvature_Quark_F2H1[2][9][0] = -((SMConstants.C_MassTop*conj(V31))/v1);
  Curvature_Quark_F2H1[2][9][1] = (II*SMConstants.C_MassTop*conj(V31))/v1;
  Curvature_Quark_F2H1[2][10][0] = -((SMConstants.C_MassTop*conj(V32))/v1);
  Curvature_Quark_F2H1[2][10][1] = (II*SMConstants.C_MassTop*conj(V32))/v1;
  Curvature_Quark_F2H1[2][11][0] = -((SMConstants.C_MassTop*conj(V33))/v1);
  Curvature_Quark_F2H1[2][11][1] = (II*SMConstants.C_MassTop*conj(V33))/v1;
  Curvature_Quark_F2H1[3][6][0] = (SMConstants.C_MassDown*V11)/v1;
  Curvature_Quark_F2H1[3][6][1] = (II*SMConstants.C_MassDown*V11)/v1;
  Curvature_Quark_F2H1[3][7][0] = (SMConstants.C_MassDown*V21)/v1;
  Curvature_Quark_F2H1[3][7][1] = (II*SMConstants.C_MassDown*V21)/v1;
  Curvature_Quark_F2H1[3][8][0] = (SMConstants.C_MassDown*V31)/v1;
  Curvature_Quark_F2H1[3][8][1] = (II*SMConstants.C_MassDown*V31)/v1;
  Curvature_Quark_F2H1[3][9][2] = (II*SMConstants.C_MassDown)/v1;
  Curvature_Quark_F2H1[3][9][3] = SMConstants.C_MassDown/v1;
  Curvature_Quark_F2H1[4][6][0] = (SMConstants.C_MassStrange*V12)/v1;
  Curvature_Quark_F2H1[4][6][1] = (II*SMConstants.C_MassStrange*V12)/v1;
  Curvature_Quark_F2H1[4][7][0] = (SMConstants.C_MassStrange*V22)/v1;
  Curvature_Quark_F2H1[4][7][1] = (II*SMConstants.C_MassStrange*V22)/v1;
  Curvature_Quark_F2H1[4][8][0] = (SMConstants.C_MassStrange*V32)/v1;
  Curvature_Quark_F2H1[4][8][1] = (II*SMConstants.C_MassStrange*V32)/v1;
  Curvature_Quark_F2H1[4][10][2] = (II*SMConstants.C_MassStrange)/v1;
  Curvature_Quark_F2H1[4][10][3] = SMConstants.C_MassStrange/v1;
  Curvature_Quark_F2H1[5][6][0] = (SMConstants.C_MassBottom*V13)/v1;
  Curvature_Quark_F2H1[5][6][1] = (II*SMConstants.C_MassBottom*V13)/v1;
  Curvature_Quark_F2H1[5][7][0] = (SMConstants.C_MassBottom*V23)/v1;
  Curvature_Quark_F2H1[5][7][1] = (II*SMConstants.C_MassBottom*V23)/v1;
  Curvature_Quark_F2H1[5][8][0] = (SMConstants.C_MassBottom*V33)/v1;
  Curvature_Quark_F2H1[5][8][1] = (II*SMConstants.C_MassBottom*V33)/v1;
  Curvature_Quark_F2H1[5][11][2] = (II*SMConstants.C_MassBottom)/v1;
  Curvature_Quark_F2H1[5][11][3] = SMConstants.C_MassBottom/v1;
  Curvature_Quark_F2H1[6][0][2] = (-II*SMConstants.C_MassUp)/v1;
  Curvature_Quark_F2H1[6][0][3] = SMConstants.C_MassUp/v1;
  Curvature_Quark_F2H1[6][3][0] = (SMConstants.C_MassDown*V11)/v1;
  Curvature_Quark_F2H1[6][3][1] = (II*SMConstants.C_MassDown*V11)/v1;
  Curvature_Quark_F2H1[6][4][0] = (SMConstants.C_MassStrange*V12)/v1;
  Curvature_Quark_F2H1[6][4][1] = (II*SMConstants.C_MassStrange*V12)/v1;
  Curvature_Quark_F2H1[6][5][0] = (SMConstants.C_MassBottom*V13)/v1;
  Curvature_Quark_F2H1[6][5][1] = (II*SMConstants.C_MassBottom*V13)/v1;
  Curvature_Quark_F2H1[7][1][2] = (-II*SMConstants.C_MassCharm)/v1;
  Curvature_Quark_F2H1[7][1][3] = SMConstants.C_MassCharm/v1;
  Curvature_Quark_F2H1[7][3][0] = (SMConstants.C_MassDown*V21)/v1;
  Curvature_Quark_F2H1[7][3][1] = (II*SMConstants.C_MassDown*V21)/v1;
  Curvature_Quark_F2H1[7][4][0] = (SMConstants.C_MassStrange*V22)/v1;
  Curvature_Quark_F2H1[7][4][1] = (II*SMConstants.C_MassStrange*V22)/v1;
  Curvature_Quark_F2H1[7][5][0] = (SMConstants.C_MassBottom*V23)/v1;
  Curvature_Quark_F2H1[7][5][1] = (II*SMConstants.C_MassBottom*V23)/v1;
  Curvature_Quark_F2H1[8][2][2] = (-II*SMConstants.C_MassTop)/v1;
  Curvature_Quark_F2H1[8][2][3] = SMConstants.C_MassTop/v1;
  Curvature_Quark_F2H1[8][3][0] = (SMConstants.C_MassDown*V31)/v1;
  Curvature_Quark_F2H1[8][3][1] = (II*SMConstants.C_MassDown*V31)/v1;
  Curvature_Quark_F2H1[8][4][0] = (SMConstants.C_MassStrange*V32)/v1;
  Curvature_Quark_F2H1[8][4][1] = (II*SMConstants.C_MassStrange*V32)/v1;
  Curvature_Quark_F2H1[8][5][0] = (SMConstants.C_MassBottom*V33)/v1;
  Curvature_Quark_F2H1[8][5][1] = (II*SMConstants.C_MassBottom*V33)/v1;
  Curvature_Quark_F2H1[9][0][0] = -((SMConstants.C_MassUp*conj(V11))/v1);
  Curvature_Quark_F2H1[9][0][1] = (II*SMConstants.C_MassUp*conj(V11))/v1;
  Curvature_Quark_F2H1[9][1][0] = -((SMConstants.C_MassCharm*conj(V21))/v1);
  Curvature_Quark_F2H1[9][1][1] = (II*SMConstants.C_MassCharm*conj(V21))/v1;
  Curvature_Quark_F2H1[9][2][0] = -((SMConstants.C_MassTop*conj(V31))/v1);
  Curvature_Quark_F2H1[9][2][1] = (II*SMConstants.C_MassTop*conj(V31))/v1;
  Curvature_Quark_F2H1[9][3][2] = (II*SMConstants.C_MassDown)/v1;
  Curvature_Quark_F2H1[9][3][3] = SMConstants.C_MassDown/v1;
  Curvature_Quark_F2H1[10][0][0] = -((SMConstants.C_MassUp*conj(V12))/v1);
  Curvature_Quark_F2H1[10][0][1] = (II*SMConstants.C_MassUp*conj(V12))/v1;
  Curvature_Quark_F2H1[10][1][0] = -((SMConstants.C_MassCharm*conj(V22))/v1);
  Curvature_Quark_F2H1[10][1][1] = (II*SMConstants.C_MassCharm*conj(V22))/v1;
  Curvature_Quark_F2H1[10][2][0] = -((SMConstants.C_MassTop*conj(V32))/v1);
  Curvature_Quark_F2H1[10][2][1] = (II*SMConstants.C_MassTop*conj(V32))/v1;
  Curvature_Quark_F2H1[10][4][2] = (II*SMConstants.C_MassStrange)/v1;
  Curvature_Quark_F2H1[10][4][3] = SMConstants.C_MassStrange/v1;
  Curvature_Quark_F2H1[11][0][0] = -((SMConstants.C_MassUp*conj(V13))/v1);
  Curvature_Quark_F2H1[11][0][1] = (II*SMConstants.C_MassUp*conj(V13))/v1;
  Curvature_Quark_F2H1[11][1][0] = -((SMConstants.C_MassCharm*conj(V23))/v1);
  Curvature_Quark_F2H1[11][1][1] = (II*SMConstants.C_MassCharm*conj(V23))/v1;
  Curvature_Quark_F2H1[11][2][0] = -((SMConstants.C_MassTop*conj(V33))/v1);
  Curvature_Quark_F2H1[11][2][1] = (II*SMConstants.C_MassTop*conj(V33))/v1;
  Curvature_Quark_F2H1[11][5][2] = (II*SMConstants.C_MassBottom)/v1;
  Curvature_Quark_F2H1[11][5][3] = SMConstants.C_MassBottom/v1;

}

bool Class_Potential_TRSM::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_TRSM::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_TRSM::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (pow(v[3],4)*lphi)/4. + pow(v[5],4)*ls + (pow(v[3],2)*(pow(v[5],2)*lphis + pow(v[4],2)*lphix + mphiSq))/2. + pow(v[5],2)*(pow(v[4],2)*lsx + msSq) + pow(v[4],2)*(pow(v[4],2)*lx + mxSq);
  
  return res;
}

double Class_Potential_TRSM::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (pow(v[3],4)*dlphi)/4. + pow(v[5],4)*dls + (pow(v[3],2)*(pow(v[5],2)*dlphis + pow(v[4],2)*dlphix + dmphiSq))/2. + pow(v[5],2)*(pow(v[4],2)*dlsx + dmsSq) + v[3]*dT4 + v[4]*(pow(v[4],3)*dlx + v[4]*dmxSq + dT5) + v[5]*dT6;
  
  return res;
}

void Class_Potential_TRSM::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
