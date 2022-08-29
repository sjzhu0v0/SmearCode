#include "TOF.h"

TOF::TOF()
{
  det_tof.push_back(new tofSmear(175, -3.7, -1.69, 20. / sqrt(2.), 20., tofSmear::endcap));
  det_tof.push_back(new tofSmear(72.7, -1.7, 1.33, 20., 20., tofSmear::barrel));
  det_tof.push_back(new tofSmear(292.1, 1.31, 4.21, 20. / sqrt(2.), 20, tofSmear::endcap));
}

int TOF::in_which_region(double eta)
{
  for (int iTof = 0; iTof < det_tof.size(); iTof++)
  {
    if (det_tof.at(iTof)->valid(eta))
    {
      return iTof + 1;
    }
  }
  return 0;
}

double TOF::get_1_beta(double p, double m, double eta)
{
  if (in_which_region(eta) == 0)
  {
    return 0.;
  }
  else
  {
    return det_tof.at(in_which_region(eta) - 1)->get_1_beta(p, m, eta);
  }
}

double TOF::get_1_beta_sigma(double eta)
{
  if (in_which_region(eta) == 0)
  {
    return 100.;
  }
  else
  {
    return det_tof.at(in_which_region(eta) - 1)->get_1_beta_sigma(eta);
  }
}