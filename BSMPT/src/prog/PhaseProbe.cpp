// Refine traced local minima at a common T; no cosmological-history inference.
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/Logger.h>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <functional>
using namespace BSMPT;
int main(int argc,char**argv) {
  if(argc!=2) { std::cerr<<"Usage: PhaseProbe point.tsv\n"; return 2; }
  Logger::SetOStream(std::cerr);
  try {
    std::ifstream in(argv[1]); std::string header,line;
    if(!std::getline(in,header)||!std::getline(in,line)) return 2;
    std::shared_ptr<Class_Potential_Origin> model=ModelID::FChoose(ModelID::ModelIDs::TRSM,GetSMConstants());
    model->setUseIndexCol(header); model->initModel(line);
    MinimumTracer tracer(model,Minimizer::WhichMinimizerDefault,false);
    std::cout<<std::setprecision(17);
    double t,a,b,c; int id;
    while(std::cin>>t>>id>>a>>b>>c) {
      std::vector<double> x{a,b,c};
      double eps=.02, scale=1+t*t;
      std::function<double(std::vector<double>)> f=[&](auto p){return model->VEff(model->MinimizeOrderVEV(p),t)/scale;};
      std::function<std::vector<double>(std::vector<double>)> df=[&](auto p){return NablaNumerical(p,f,eps);};
      std::function<std::vector<std::vector<double>>(std::vector<double>)> hf=[&](auto p){return HessianNumerical(p,f,eps);};
      try {
        x=tracer.LocateMinimum(x,df,hf,1e-7,1e-2,100);
        double v_coarse=f(x)*scale; eps=.01;
        x=tracer.LocateMinimum(x,df,hf,1e-7,1e-2,100);
        for(auto &v:x) v=std::abs(v); // gauge / exact discrete symmetry copies
        double v=f(x)*scale, grad=0;
        auto g=df(x); for(double q:g) grad+=q*q;
        grad=std::sqrt(grad); auto h=hf(x);
        Eigen::Matrix3d H; for(int i=0;i<3;i++)for(int j=0;j<3;j++)H(i,j)=h[i][j]*scale;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(H);
        double eigen=solver.eigenvalues().minCoeff();
        double error=std::max({1e-5,std::abs(v)*2e-12,4*std::abs(v-v_coarse)});
        bool good=std::isfinite(v)&&std::isfinite(grad)&&grad<2e-5&&eigen>1e-5;
        std::cout<<"PROBE "<<t<<' '<<id<<' '<<(good?"minimum":(std::isfinite(grad)&&grad<2e-5&&eigen< -1e-5?"saddle":"unresolved"))<<' '<<x[0]<<' '<<x[1]<<' '<<x[2]<<' '<<v<<' '<<error<<' '<<grad*scale<<' '<<eigen<<std::endl;
      } catch(const std::exception &e) { std::cout<<"PROBE "<<t<<' '<<id<<" error 0 0 0 0 0 0 0"<<std::endl; std::cerr<<e.what()<<'\n'; }
    }
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
