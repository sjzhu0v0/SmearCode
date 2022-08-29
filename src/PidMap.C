#include "PidMap.h"

void PidMap::GetPidMap(PID_SYSTEM pid_system, double *pidmap) {
  int PidIndex = 0;
  for (int iParticleType = 0; iParticleType < n_particle_type;
       iParticleType++) {
    for (int iDetType = 0; iDetType < n_particle_type; iDetType++) {
      *(pidmap + iParticleType*n_det_type + iDetType) = 10;
    }
  }

  if (pid_system == EICC_ACCEPTANCE) {
    int region_index = 0;
    for (int iRegion = 0; iRegion < 3; iRegion++) {
      if (pseduorapidity > min_pid_eta_region[iRegion] &&
          pseduorapidity < max_pid_eta_region[iRegion]) {
        region_index = iRegion;
      }
    }
    if (momentum < PID_acceptance[pdg_id][region_index]) {
      PidIndex = (int)pdg_id;
      for (int iDetType = 0; iDetType < n_det_type; iDetType++) {
        *(pidmap + PidIndex*n_det_type + iDetType) = 0;
      }
    } else {
      for (int iParticleType = 0; iParticleType < n_particle_type;
           iParticleType++) {
        for (int iDetType = 0; iDetType < n_particle_type; iDetType++) {
          *(pidmap + iParticleType*n_det_type + iDetType) = 1.;
        }
      }
    }
  } else {
    for (int iParticleType = 0; iParticleType < n_particle_type;
         iParticleType++) {
      for (int iDetType = 0; iDetType < n_particle_type; iDetType++) {
        *(pidmap + iParticleType*n_det_type + iDetType) =
            Nsigma_to_particle[iParticleType][iDetType];
      }
    }
  }
}

int PidMap::get_expected_pid_int(PID_SYSTEM pid_system) {
  int PidIndex = 0;
  if (pid_system == EICC_ACCEPTANCE) {
    int region_index = 0;
    for (int iRegion = 0; iRegion < 3; iRegion++) {
      if (pseduorapidity > min_pid_eta_region[iRegion] &&
          pseduorapidity < max_pid_eta_region[iRegion]) {
        region_index = iRegion;
      }
    }
    if (momentum < PID_acceptance[pdg_id][region_index]) {
      PidIndex = (int)pdg_id;
    } else {
      return 1;
      // returned value 1 means that the chosen PID system has no
      // identification power with the input mc info.
    }
  } else {
    cal_pidmap();
    if (prob_sum == 0) {
      return 1;
    }
    for (int i = 0; i < n_particle_type; i++) {
      if (bayes_probability_to_particle[i] >
          bayes_probability_to_particle[PidIndex]) {
        PidIndex = i;
      }
    }
  }
  return PDGID[PidIndex];
}

PidMap::PidMap(PidMap::PID pid, double eta, double p) {
  pdg_id = pid;
  momentum = p;
  if (p > 20) {
    p = 19.9;
  } else {
    momentum = p;
  }
  pseduorapidity = eta;
  if (eta < -3.5) {
    pseduorapidity = -3.499;
  }
  if (eta > 3.5) {
    pseduorapidity = 3.499;
  }
  smear_particle_mcinfo();
}

PidMap::PidMap(int pid, double eta, double p) {
  for (int i = 0; i < n_particle_type; i++) {
    if (PDGID[i] == abs(pid)) {
      pdg_id = all_particle_type[i];
      break;
    }
  }
  momentum = p;
  if (p > 20) {
    p = 19.9;
  } else {
    momentum = p;
  }
  pseduorapidity = eta;
  if (eta < -3.5) {
    pseduorapidity = -3.499;
  }
  if (eta > 3.5) {
    pseduorapidity = 3.499;
  }
  smear_particle_mcinfo();
}

void PidMap::is_hit() {
  for (int i = 0; i < n_det_type; i++) {
    det_hit[i] = in_which_region(i);
  }
  // cherenkov detector: checking whether threshold is met

  if (momentum > min_momentum_applying_emc[det_hit[2] - 1] && det_hit[2] != 0) {
    isEmcUsed = true;
  }

  if (det_hit[1] == 1 && momentum < momentum_threshold_mRICH(pdg_id))
    det_hit[1] = -1;
  if (det_hit[1] == 2 && momentum < momentum_threshold_hpDIRC(pdg_id))
    det_hit[1] = -2;
  if (det_hit[1] == 3 && momentum < momentum_threshold_dRICH(0, pdg_id))
    det_hit[1] = -3;
}

int PidMap::in_which_region(int iDet) {
  if (iDet == 0) {
    for (int i = 0; i < 3; i++) {
      if (pseduorapidity > min_tof_eta_region[i] &&
          pseduorapidity < max_tof_eta_region[i]) {
        return i + 1;
      }
    }
  } else if (iDet == 1) {
    for (int i = 0; i < 3; i++) {
      if (pseduorapidity > min_cherenkov_eta_region[i] &&
          pseduorapidity < max_cherenkov_eta_region[i]) {
        return i + 1;
      }
    }
  } else if (iDet == 2) {
    for (int i = 0; i < 3; i++) {
      if (pseduorapidity > min_emc_eta_region[i] &&
          pseduorapidity < max_emc_eta_region[i]) {
        return i + 1;
      }
    }
    return 0;
  } else {
    cerr << "iDet is wrong" << endl;
    return -10;
  }
  return 0;
}

void PidMap::cal_pidmap() {
  // ECCE PID system is used to identify the input particles.
  for (int iPID = 0; iPID < n_particle_type; iPID++) {
    for (int iDet = 0; iDet < n_det_type; iDet++) {
      TotalSigma_to_particle[iPID] += pow(Nsigma_to_particle[iPID][iDet], 2);
    }
    TotalSigma_to_particle[iPID] = sqrt(TotalSigma_to_particle[iPID]);
  }

  for (int iPID = 0; iPID < n_particle_type; iPID++) {
    for (int iDet = 0; iDet < n_det_type; iDet++) {
      if (!isEmcUsed && det_hit[1] < 1) {
        bayes_probability_to_particle[iPID] *= ROOT::Math::chisquared_cdf_c(
            pow(TotalSigma_to_particle[iPID], 2), 1.);
      } else if (!isEmcUsed && det_hit[1] > 1) {
        bayes_probability_to_particle[iPID] *= ROOT::Math::chisquared_cdf_c(
            pow(TotalSigma_to_particle[iPID], 2), 2.);
      } else {
        bayes_probability_to_particle[iPID] *= ROOT::Math::chisquared_cdf_c(
            pow(TotalSigma_to_particle[iPID], 2), 3.);
      }
    }
  }

  prob_sum = 0;

  for (int iPID = 0; iPID < n_particle_type; iPID++)
    prob_sum += bayes_probability_to_particle[iPID];

  if (prob_sum == 0) {
    cout << "prob_sum is zero. PidMap info: " << fMass[pdg_id] << " "
         << momentum << " " << pseduorapidity << endl;
    for (int iPID = 0; iPID < n_particle_type; iPID++)
      bayes_probability_to_particle[iPID] /= 1.;
  } else {
    for (int iPID = 0; iPID < n_particle_type; iPID++)
      bayes_probability_to_particle[iPID] /= prob_sum;
  }
}

void PidMap::smear_particle_mcinfo() {
  is_hit();
  for (int iPID = 0; iPID < n_particle_type; iPID++) {
    Nsigma_to_particle[iPID][0] =
        Nsigma_to_particle_tof(all_particle_type[iPID]);
    Nsigma_to_particle[iPID][1] =
        Nsigma_to_particle_cherenkov(all_particle_type[iPID]);
    Nsigma_to_particle[iPID][2] =
        Nsigma_to_particle_emc(all_particle_type[iPID]);
  }
}

double PidMap::probability_to_particle_tof(PidMap::PID toWhichParticle) {
  return 1 - 2 * abs(ROOT::Math::normal_cdf_c(
                         Nsigma_to_particle_cherenkov(toWhichParticle)) -
                     0.5);
}

double PidMap::sigma_to_particle_tof() {
  return MTof.get_1_beta_sigma(pseduorapidity);
}

double PidMap::Nsigma_to_particle_tof(PidMap::PID toWhichParticle) {
  double tofTrue_1Beta =
      MTof.get_1_beta(momentum, fMass[pdg_id], pseduorapidity);
  double tofExpected_1Beta =
      MTof.get_1_beta(momentum, fMass[toWhichParticle], pseduorapidity);
  double tof_1BetaSigma = MTof.get_1_beta_sigma(pseduorapidity);
  return (tofTrue_1Beta + rnd_tof * tof_1BetaSigma - tofExpected_1Beta) /
         tof_1BetaSigma;
}

double PidMap::probability_to_particle_cherenkov(PidMap::PID toWhichParticle) {
  if (det_hit[1] == 0) {
    return 0.5;
  }
  if (det_hit[1] == 1) {
    if (momentum > momentum_threshold_mRICH(toWhichParticle)) {
      return probability_to_particle_mRICH(toWhichParticle);
    } else {
      return 0.;
    }
  }
  if (det_hit[1] == 2) {
    if (momentum > momentum_threshold_hpDIRC(toWhichParticle)) {
      return probability_to_particle_hpDIRC(toWhichParticle);
    } else {
      return 0.;
    }
  }
  if (det_hit[1] == 3) {
    if (momentum > momentum_threshold_dRICH(0, toWhichParticle)) {
      if (momentum > momentum_threshold_dRICH(1, toWhichParticle)) {
        return probability_to_particle_dRICH_C2F6(toWhichParticle);
      } else {
        return probability_to_particle_dRICH_aerogel(toWhichParticle);
      }
    } else {
      return 0.;
    }
  }
  if (det_hit[1] == -1) {
    if (momentum < momentum_threshold_dRICH(0, toWhichParticle)) {
      return 0.5;
    } else {
      return 0.;
    }
  }
  if (det_hit[1] == -2) {
    if (momentum < momentum_threshold_hpDIRC(toWhichParticle)) {
      return 0.5;
    } else {
      return 0.;
    }
  }
  if (det_hit[1] == -3) {
    if (momentum < momentum_threshold_mRICH(toWhichParticle)) {
      return 0.5;
    } else {
      return 0.;
    }
  }
  return -1.;
}

double PidMap::Nsigma_to_particle_cherenkov(PidMap::PID toWhichParticle) {
  if (det_hit[1] == 0) {
    return 0.;
  } else if (det_hit[1] == 1) {
    if (momentum > momentum_threshold_mRICH(toWhichParticle)) {
      return Nsigma_to_particle_mRICH(toWhichParticle);
    } else {
      return 100.;
    }
  } else if (det_hit[1] == 2) {
    if (momentum > momentum_threshold_hpDIRC(toWhichParticle)) {
      return Nsigma_to_particle_hpDIRC(toWhichParticle);
    } else {
      return 100.;
    }
  } else if (det_hit[1] == 3) {
    if (momentum > momentum_threshold_dRICH(0, toWhichParticle)) {
      if (momentum > momentum_threshold_dRICH(1, pdg_id)) {
        if (momentum > momentum_threshold_dRICH(1, toWhichParticle)) {
          return Nsigma_to_particle_dRICH_C2F6(toWhichParticle);
        } else {
          return 100;
        }
      } else {
        return Nsigma_to_particle_dRICH_aerogel(toWhichParticle);
      }
    } else {
      return 100.;
    }
  } else if (det_hit[1] == -1) {
    if (momentum < momentum_threshold_mRICH(toWhichParticle)) {
      return 0.;
    } else {
      return 100.;
    }
  } else if (det_hit[1] == -2) {
    if (momentum < momentum_threshold_hpDIRC(toWhichParticle)) {
      return 0.;
    } else {
      return 100.;
    }
  } else if (det_hit[1] == -3) {
    if (momentum < momentum_threshold_dRICH(0, toWhichParticle)) {
      return 0.;
    } else {
      return 100.;
    }
  }
  return -1.;
}

double PidMap::probability_to_particle_emc(PidMap::PID toWhichParticle) {
  if (det_hit[2] == 0) {
    return 1.;
  }

  if (momentum < min_momentum_applying_emc[det_hit[2] - 1]) {
    return 0.5;
  } else {
    if (pdg_id == PidMap::electron) {
      if (toWhichParticle == PidMap::electron) {
        if (rnd_emc < 0.95) {
          return 1.;
        } else {
          return 0.;
        }
      } else {
        if (rnd_emc < 0.95) {
          return 0.;
        } else {
          return 1.;
        }
      }
    } else if (pdg_id == PidMap::pion) {
      double pi_survive =
          1. / pow(10, gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum)));
      if (toWhichParticle == PidMap::electron) {
        if (rnd_emc < pi_survive) {
          return 1.;
        } else {
          return 0.;
        }
      } else {
        if (rnd_emc < pi_survive) {
          return 0.;
        } else {
          return 1.;
        }
      }
    }
  }
  return -1.;
}

double PidMap::Nsigma_to_particle_emc(PidMap::PID toWhichParticle) {
  if (det_hit[2] == 0) {
    return 0.;
  } else {
    if (momentum < min_momentum_applying_emc[det_hit[2] - 1]) {
      return 0.;
    } else {
      if (pdg_id == PidMap::electron) {
        if (toWhichParticle == PidMap::electron) {
          if (rnd_emc < 0.95) {
            return 0.;
          } else {
            return gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum));
          }
        } else {
          if (rnd_emc < 0.95) {
            return gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum));
          } else {
            return 0.;
          }
        }
      } else {
        double pi_survive =
            1. /
            pow(10, gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum)));
        if (toWhichParticle == PidMap::electron) {
          if (rnd_emc < pi_survive) {
            return 0.;
          } else {
            return gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum));
          }
        } else {
          if (rnd_emc < pi_survive) {
            return gr_emc_log_pi_rej[det_hit[2] - 1]->Eval(log(momentum));
          } else {
            return 0.;
          }
        }
      }
    }
  }
}

double PidMap::probability_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle) {
  return 1 - 2 * abs(ROOT::Math::normal_cdf_c(
                         Nsigma_to_particle_dRICH_C2F6(toWhichParticle)) -
                     0.5);
}

double PidMap::sigma_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle) {
  return Mdual_rich_c2f6.cherenkovAngleSigma(eta_to_theta(pseduorapidity),
                                             momentum, fMass[toWhichParticle]);
}

double PidMap::Nsigma_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle) {
  double temp_debug =
      sqrt(fMass[pdg_id] * fMass[pdg_id] + momentum * momentum) /
      (1.0014 * momentum);
  double cherenkovTrueAngle =
      Mdual_rich_c2f6.cherenkovAngle(momentum, fMass[pdg_id]);
  double cherenkovExpectedAngle =
      Mdual_rich_c2f6.cherenkovAngle(momentum, fMass[toWhichParticle]);
  return (cherenkovTrueAngle - cherenkovExpectedAngle +
          rnd_cherenkov * sigma_to_particle_dRICH_C2F6(pdg_id)) /
         sigma_to_particle_dRICH_C2F6(toWhichParticle);
}

double
PidMap::probability_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle) {
  return 1 - 2 * abs(ROOT::Math::normal_cdf_c(
                         Nsigma_to_particle_dRICH_aerogel(toWhichParticle)) -
                     0.5);
}

double PidMap::sigma_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle) {
  return Mdual_rich_aerogel.cherenkovAngleSigma(
      eta_to_theta(pseduorapidity), momentum, fMass[toWhichParticle]);
}

double PidMap::Nsigma_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle) {
  double cherenkovTrueAngle =
      Mdual_rich_aerogel.cherenkovAngle(momentum, fMass[pdg_id]);
  double cherenkovExpectedAngle =
      Mdual_rich_aerogel.cherenkovAngle(momentum, fMass[toWhichParticle]);
  return (cherenkovTrueAngle - cherenkovExpectedAngle +
          rnd_cherenkov * sigma_to_particle_dRICH_aerogel(pdg_id)) /
         sigma_to_particle_dRICH_aerogel(toWhichParticle);
}

double PidMap::probability_to_particle_hpDIRC(PidMap::PID toWhichParticle) {
  return 1 - 2 * abs(ROOT::Math::normal_cdf_c(
                         Nsigma_to_particle_hpDIRC(toWhichParticle)) -
                     0.5);
}

double PidMap::sigma_to_particle_hpDIRC(PidMap::PID toWhichParticle) {
  return MdrcPidFast.GetCctr(0, eta_to_theta(pseduorapidity), momentum, 0.5);
}

double PidMap::Nsigma_to_particle_hpDIRC(PidMap::PID toWhichParticle) {
  double cherenkovTrueAngle = MdrcPidFast.GetAng(fMass[pdg_id], momentum);
  double cherenkovExpectedAngle =
      MdrcPidFast.GetAng(fMass[toWhichParticle], momentum);
  return (cherenkovTrueAngle - cherenkovExpectedAngle +
          rnd_cherenkov * sigma_to_particle_hpDIRC(pdg_id)) /
         sigma_to_particle_hpDIRC(toWhichParticle);
}

double PidMap::probability_to_particle_mRICH(PidMap::PID toWhichParticle) {
  return 1 - 2 * abs(ROOT::Math::normal_cdf_c(
                         Nsigma_to_particle_mRICH(toWhichParticle)) -
                     0.5);
}

double PidMap::sigma_to_particle_mRICH(PidMap::PID toWhichParticle) {
  return MmRICH.getdAng(fMass[toWhichParticle], momentum);
}

double PidMap::Nsigma_to_particle_mRICH(PidMap::PID toWhichParticle) {
  double cherenkovTrueAngle = MmRICH.getAng(fMass[pdg_id], momentum);
  double cherenkovExpectedAngle =
      MmRICH.getAng(fMass[toWhichParticle], momentum);
  // cout << "Nsigma_to_particle_mRICH: " << cherenkovTrueAngle << " "
  //      << cherenkovExpectedAngle << " " << rnd_cherenkov << " "
  //      << sigma_to_particle_mRICH(toWhichParticle) << endl;
  return (cherenkovTrueAngle - cherenkovExpectedAngle +
          rnd_cherenkov * sigma_to_particle_mRICH(pdg_id)) /
         sigma_to_particle_mRICH(toWhichParticle);
}

// type1: threshold of c2f6 is used. type0: threshold of aerogel is used.
double PidMap::momentum_threshold_dRICH(int type, PidMap::PID toWhichParticle) {
  if (type == 0) {
    return Mdual_rich_aerogel.cherenkovThreshold(fMass[toWhichParticle]);
  } else if (type == 1) {
    return Mdual_rich_c2f6.cherenkovThreshold(fMass[toWhichParticle]);
  } else {
    return 100.;
  }
}

double PidMap::momentum_threshold_hpDIRC(PidMap::PID toWhichParticle) {
  return MdrcPidFast.MomentumThreshold(fMass[toWhichParticle]);
}

double PidMap::momentum_threshold_mRICH(PidMap::PID toWhichParticle) {
  return MmRICH.MomentumThreshold(fMass[toWhichParticle]);
}
