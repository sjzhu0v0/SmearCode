#include "tofSmear.h"

#include "string"
#include "vector"

#ifndef __TOF_H__
#define __TOF_H__

class TOF
{
private:
    vector<tofSmear *> det_tof;

public:
    TOF();
    virtual ~TOF() {}
    
    int in_which_region(double eta);
    double get_1_beta(double p, double m, double eta);
    double get_1_beta_sigma(double eta);
};

#endif
