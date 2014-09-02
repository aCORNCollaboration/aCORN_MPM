#ifndef BASEDATASCANNER_HH
#define BASEDATASCANNER_HH

#include "Enums.hh"
#include <Rtypes.h>
#include "RunSetScanner.hh"

const unsigned int NCH_MAX = 32;        ///< number of channels in DAQ
const unsigned int N_E_PMT = 19;        ///< number of main electron scintillator PMTs
const unsigned int N_V_PMT = 8;         ///< number of veto PMTs
const unsigned int NS_PER_CLOCK = 10;   ///< detector time-of-flight clock pulse

/// Generic class for processed data TChains
class BaseDataScanner: public RunSetScanner {
public:
    /// Constructor
    BaseDataScanner(const std::string& treeName, bool fp);
    
    bool is4p;                  ///< whether this is "4-proton" style data
    
    Long_t T_p;                 ///< proton (energy pulse) arrival time, loaded in 10ns units, calibrated to [ns]
    Long_t T_d;                 ///< [4p] discriminator arrival time
    Int_t nPSig;                ///< [4p] number of proton signals
    Int_t T_e2p;                ///< proton time-of-flight from electron signal, loaded in 10ns, calibrated to [ns]
    Int_t T_e2d;                ///< [4p] proton discriminator time to electron signal
    
    Int_t E_p;                  ///< proton energy (uncalibrated channels)
    Int_t E_p_0;                ///< previous proton energy for multiple-electron flag
    Int_t E_e;                  ///< Electron energy (uncalibrated channels)
    Double_t E_recon;           ///< Calibrated reconstructed electron energy
    
    Int_t nP;                   ///< number of protons between electron trigger and this proton event
    UInt_t nE;                  ///< number of electron PMTs triggered
    Int_t V;                    ///< number of veto PMTs triggered
    
    UInt_t DetFired;            ///< detector firing bit map 0-30
    UInt_t DetPiled;            ///< pileup flags [e- dead zone][det. 0-30]
    Int_t Idx;                  ///< proton index
    Short_t E_PMT[NCH_MAX];     ///< individual PMT energy
    Char_t T_PMT[NCH_MAX];      ///< individual PMT time
    Char_t T_PMT_median;        ///< median PMT timing offset
    Char_t Max_PMT;             ///< PMT with largest signal
    Short_t E_Max_PMT;          ///< energy of max PMT signal
    
    TriggerCategory flags;      ///< Event category flags
    
    /// generate event flags
    virtual void makeFlags();
    /// load calibrations for new run
    virtual void loadNewRun(RunID r);
    /// calibrations after loading event
    virtual void calibrate();
    
    
    double physicsWeight;       ///< simulated event spectrum re-weighting factor
   
protected:
    std::map<RunID, AcornCalibrator*> cals;     ///< calibrators for each run
    AcornCalibrator* currentCal;                ///< calibrator for current run
    
    /// load calibrator for current run
    void loadCal(RunID rn);
};

#endif
