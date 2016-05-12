#ifndef GENANALYZER_H
#define GENANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TH3.h"
#include "TKey.h"

class GenAnalyzer {
    public:
        GenAnalyzer();
        ~GenAnalyzer();
        //virtual float GetDYWeight(const edm::Event&);
        virtual float GetPUWeight(const edm::Event&);
    //    virtual float GetPDFWeight(const edm::Event&);
        virtual std::pair<float, float> GetQ2Weight(const edm::Event&);

      
    private:
        /*
        bool isDYFile;
        int npMax;
        float ptMax, htMax;
        std::string Sample;
        
        TH1* Sum;
        TH1* Num;
        TFile* DYFile;
        */
        edm::LumiReWeighting* LumiWeights;
};


//namespace LHAPDF {
//  void initPDFSet(int nset, int setid, int member=0);
//  void initPDFSet(int nset, const std::string& filename, int member=0);
//  int numberPDF(int nset);
//  void usePDFMember(int nset, int member);
//  double xfx(int nset, double x, double Q, int fl);
//  double getXmin(int nset, int member);
//  double getXmax(int nset, int member);
//  double getQ2min(int nset, int member);
//  double getQ2max(int nset, int member);
//  void extrapolate(bool extrapolate=true);
//}


#endif
