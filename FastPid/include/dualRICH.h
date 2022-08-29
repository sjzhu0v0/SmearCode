/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   April 2020

#ifndef __DUALRICH_H__
#define __DUALRICH_H__

#include "genericRICH.h"

/** dualRICH aerogel **/

class dualRICH_aerogel : public genericRICH
{
 public:
  dualRICH_aerogel();
  virtual ~dualRICH_aerogel() = default;
};


/** dualRICH C2F6 **/

class dualRICH_C2F6 : public genericRICH
{
 public:
  dualRICH_C2F6();
  virtual ~dualRICH_C2F6() = default;
};

#endif /* __DUALRICH_H__ */