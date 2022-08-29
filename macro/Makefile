PROGRAM  = Particle_Mc_info
ROOTLIBS = $(shell root-config --libs --cflags)
CXX      =  g++
FASTPID = ../FastPid
PACKAGE = ../
EICCSMEAR_LIB = /home/szhu/Software/EiccSoft/EicC_Mvd_DP/src/libEicC_Mvd_DP.so
PYTHIA8_LIB = /home/szhu/Software/PYTHIA8/lib/libpythia8.so
EICSMEAR_LIB = /home/szhu/Software/EicSoft/lib/libeicsmear.so 
EXTERNAL_LIB = $(EICSMEAR_LIB) $(EICCSMEAR_LIB) $(shell root-config --libs) $(PYTHIA8_LIB)
CPPFLAGS =  -I$(PYTHIAPATH)/include -I$(ROOTSYS)/include -I$(LHAPDFPATH)/include -I$(PACKAGE)/include -I$(FASTPID)/include $(shell root-config --cflags) 
LDFLAGS  =  $(ROOTLIBS) 
# LDFLAGS  = -L$(PYTHIA_LIB) -Wl,-rpath,$(PYTHIA_LIB) -lpythia8 -ldl -lLHAPDF  $(ROOTLIBS) 

pid_power: pid_power.cpp PidMap.o
	$(CXX) -o ../bin/pid_power pid_power.cpp PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB)

all: cut_cal_d0 ParticleDecay Particle_Mc_info smear

cut_cal_d0_debug: cut_cal_d0_debug.cpp SmearEvent.o
	$(CXX) -o ../bin/cut_cal_d0_debug cut_cal_d0_debug.cpp SmearEvent.o PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB)

cut_cal_lambda_c: cut_cal_lambda_c.cpp SmearEvent.o
	$(CXX) -o ../bin/cut_cal_lambda_c cut_cal_lambda_c.cpp SmearEvent.o PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB)

cut_cal_d0: cut_cal_d0.cpp SmearEvent.o
	$(CXX) -o ../bin/cut_cal_d0 cut_cal_d0.cpp SmearEvent.o PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB)

ParticleDecay: McParticles_byEvent.o
	$(CXX) -o ../bin/ParticleDecay ParticleDecay.cpp SmearEvent.o PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o McParticles_byEvent.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB) $(shell pythia8-config --cflags --libs)

Particle_Mc_info: Particle_Mc_info.cpp
	$(CXX) -o ../bin/Particle_Mc_info Particle_Mc_info.cpp $(EICSMEAR_LIB) $(LDFLAGS) $(CPPFLAGS) 

Hadron_Mc_info: Hadron_Mc_info.cpp
	$(CXX) -o ../bin/Hadron_Mc_info Hadron_Mc_info.cpp $(EICSMEAR_LIB) $(LDFLAGS) $(CPPFLAGS) 

smear: smear.cpp SmearEvent.o
	$(CXX) -o ../bin/smear smear.cpp SmearEvent.o PidMap.o dualRICH.o DrcPidFast.o mRICH.o TOF.o tofSmear.o MDefinition.o genericRICH.o genericDetector.o $(CPPFLAGS) $(LDFLAGS) $(EXTERNAL_LIB)

McParticles_byEvent.o: SmearEvent.o ../src/McParticles_byEvent.C ../include/McParticles_byEvent.h
	$(CXX) -c ../src/McParticles_byEvent.C $(CPPFLAGS) $(shell pythia8-config --cflags)

SmearEvent.o: PidMap.o ../src/SmearEvent.C
	$(CXX) -c ../src/SmearEvent.C $(CPPFLAGS) 

PidMap.o: dualRICH.o DrcPidFast.o mRICH.o TOF.o MDefinition.o
	$(CXX) -c ../src/PidMap.C $(CPPFLAGS)

MDefinition.o: dualRICH.o DrcPidFast.o mRICH.o TOF.o
	$(CXX) -c ../src/MDefinition.C $(CPPFLAGS) 

dualRICH.o: genericRICH.o
	$(CXX) -c $(FASTPID)/src/dualRICH.C $(CPPFLAGS) 

genericRICH.o: genericDetector.o
	$(CXX) -c $(FASTPID)/src/genericRICH.C $(CPPFLAGS)

genericDetector.o:
	$(CXX) -c $(FASTPID)/src/genericDetector.C $(CPPFLAGS)

DrcPidFast.o: 
	$(CXX) -c $(FASTPID)/src/DrcPidFast.cxx $(CPPFLAGS) 

mRICH.o:
	$(CXX) -c $(FASTPID)/src/mRICH.C $(CPPFLAGS) 

TOF.o: tofSmear.o
	$(CXX) -c $(FASTPID)/src/TOF.C $(CPPFLAGS) 

tofSmear.o:
	$(CXX) -c $(FASTPID)/src/tofSmear.C $(CPPFLAGS) 
