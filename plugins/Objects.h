#ifndef OBJECTS_H
#define OBJECTS_H

struct LeptonType {
    LeptonType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), charge(0), pdgId(0), pfIso03(-1.), pfIso04(-1.), trkIso(-1.), miniIso(-1.), dxy(-99.), dz(-99.), dPhi_met(-1.), isElectron(false), isMuon(false), isVeto(false), isLoose(false), isMedium(false), isTight(false), isHighpt(false), isMatched(false) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    int charge;
    int pdgId;
    float pfIso03;
    float pfIso04;
    float trkIso;
    float miniIso;
    float dxy;
    float dz;
    float dPhi_met;
    bool isElectron;
    bool isMuon;
    bool isVeto;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isHighpt;
    bool isMatched;
};

struct JetType {
  JetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), ptRaw(-1.), ptUnc(-1.), dPhi_met(-1.), dPhi_jet1(-1.), puId(-1.), CSV(-99.), CSVR(-99.), chf(-1.), nhf(-1.), phf(-1.), elf(-1.), muf(-1.), chm(-1), npr(-1), flavour(0), mother(0), isLoose(false), isMedium(false), isTight(false), isCSVL(false), isCSVM(false), isCSVT(false), isMatched(false) {}
  float pt;
  float eta;
  float phi;
  float mass;
  float energy;
  float ptRaw;
  float ptUnc;
  float dPhi_met;
  float dPhi_jet1;
  float puId;
  float CSV;
  float CSVR;
  float chf;
  float nhf;
  float phf;
  float elf;
  float muf;
  int chm;
  int npr;
  int flavour;
  int mother;
  bool isLoose;
  bool isMedium;
  bool isTight;
  bool isCSVL;
  bool isCSVM;
  bool isCSVT;
  bool isMatched;
};

struct CandidateType {
  CandidateType(): pt(-1.), eta(-9.), phi(-9.), et(-1.), p(-1.), energy(-1.), ptGen(-1.), etaGen(-9.), phiGen(-9.), mass(-1.), massGen(-1.), massGenJet(-1.), dR(-1.), dEta(-1.), dPhi(-1.), twist(-1.), isMatched(false) {} // angle(-1.), ptBalance(-9.), centrality(-1.), charge(0), 
  float pt;
  float eta;
  float phi;
  float et;
  float p;
  float energy;
  float ptGen;
  float etaGen;
  float phiGen;
  float mass;
  float massGen;
  float massGenJet;
  float dR;
  float dEta;
  float dPhi;
  float twist;
//  float angle;
//  float ptBalance;
//  float centrality;
//  int charge;
  bool isMatched;
};

struct MEtType {
  MEtType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.) {}
  float pt;
  float eta;
  float phi;
  float sign;
};

struct MEtCorType {
  MEtCorType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptScaleUp(-1.), ptScaleDown(-1.), ptResUp(-1.), ptResDown(-1.) {}
  float pt;
  float eta;
  float phi;
  float sign;
  float ptScaleUp;
  float ptScaleDown;
  float ptResUp;
  float ptResDown;
};

struct MEtFullType {
  MEtFullType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptRaw(-1.), phiRaw(-1.), ptGen(-1.), phiGen(-1.), ptJERUp(-1.), ptJERDown(-1.), ptJESUp(-1.), ptJESDown(-1.), ptMUSUp(-1.), ptMUSDown(-1.), ptELSUp(-1.), ptELSDown(-1.), ptTAUUp(-1.), ptTAUDown(-1.), ptUNCUp(-1.), ptUNCDown(-1.), ptPHOUp(-1.), ptPHODown(-1.), phf(-1.), nhf(-1.), elf(-1.), chf(-1.), muf(-1.) {}
  float pt;
  float eta;
  float phi;
  float sign;
  float ptRaw;
  float phiRaw;
  float ptGen;
  float phiGen;
  float ptJERUp;
  float ptJERDown;
  float ptJESUp;
  float ptJESDown;
  float ptMUSUp;
  float ptMUSDown;
  float ptELSUp;
  float ptELSDown;
  float ptTAUUp;
  float ptTAUDown;
  float ptUNCUp;
  float ptUNCDown;
  float ptPHOUp;
  float ptPHODown;
  float phf;
  float nhf;
  float elf;
  float chf;
  float muf;
};



struct LorentzType {
  LorentzType(): pt(-1.), eta(-9.), phi(-9.), energy(-1.), mass(-1.) {}
  float pt;
  float eta;
  float phi;
  float energy;
  float mass;
};

struct EventType {
  EventType(): id1(-9), id2(-9), x1(-1.), x2(-1.), Q(-1.) {}
  int id1;
  int id2;
  float x1;
  float x2;
  float Q;
};

#endif

