#ifndef SRC_CLASSPOTENTIALTRSM_H_
#define SRC_CLASSPOTENTIALTRSM_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_Potential_TRSM : public Class_Potential_Origin
{
public:
  Class_Potential_TRSM(const ISMConstants &smConstants);
  virtual ~Class_Potential_TRSM();

  // Initialize input parameters
  double m1 = 0;
  double m2 = 0;
  double m3 = 0;
  double vs = 0;
  double a12 = 0;
  double lx = 0;
  double lphix = 0;
  double lsx = 0;

  // Initialize dependent parameters
  double v1 = 0;
  double lphi = 0;
  double ls = 0;
  double lphis = 0;
  double mxSq = 0;
  double mphiSq = 0;
  double msSq = 0;

  // Initialize counter terms
  double dmphiSq = 0;
  double dlphi = 0;
  double dmsSq = 0;
  double dls = 0;
  double dmxSq = 0;
  double dlx = 0;
  double dlphis = 0;
  double dlphix = 0;
  double dlsx = 0;
  double dT1 = 0;
  double dT2 = 0;
  double dT3 = 0;
  double dT4 = 0;
  double dT5 = 0;
  double dT6 = 0;


  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void AdjustRotationMatrix() override;
  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_TRSM_H_ */
