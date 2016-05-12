#include "BTagInterface.h"

  
BTagInterface::BTagInterface(std::string bTag) {
  Tagger=bTag;
  isBDFile=isShapesFile=false;
  std::string tag("");
  if(Tagger=="combinedSecondaryVertexBJetTags") tag="CSV";
  else if(Tagger=="combinedSecondaryVertexV1BJetTags") tag="CSVV1";
  else if(Tagger=="combinedSecondaryVertexSoftPFLeptonV1BJetTags") tag="CSVSL";
  else {
    std::cout << " - BTagInterface Error: unrecognized tagger " << tag << std::endl;
    return;
  }
  std::string wp[3]={"L", "M", "T"};
  std::string fl[3]={"B", "C", "L"};
  
  BDFile=new TFile("data/BTagPerformanceDBWinter13.root", "READ");
  BDFile->cd();
  if(!BDFile->IsZombie()) {
  	for(int i=0; i<3; i++) {
  		SFb[i]=(TH2F*)BDFile->Get(("BTagSF/"+tag+wp[i]+"_SFb").c_str());
  		SFl[i]=(TH2F*)BDFile->Get(("BTagSF/"+tag+wp[i]+"_SFl").c_str());
  		isBDFile=true;
  	}
  }
  else std::cout << " - BTagInterface Warning: No BTag Scale Factor File" << std::endl;
  
  ShapesFile=new TFile("data/BTagShapes.root", "READ");
  ShapesFile->cd();
  if(!ShapesFile->IsZombie()) {
  	for(int f=0; f<3; f++) {
  		Shape[f]=(TH1F*)ShapesFile->Get(("j_bTagDiscr"+fl[f]).c_str());
  		//Shape[f]->Smooth(100); // Smoothing
  		//Shape[f]->Rebin(10); // Rebinning
  		isShapesFile=true;
  	}
  }
  else std::cout << " - BTagInterface Warning: No Shape File" << std::endl;
  
  
  WorkingPoint[0]=0.;
  if(Tagger=="combinedSecondaryVertexBJetTags") {
    WorkingPoint[1]=0.244;
    WorkingPoint[2]=0.679;
    WorkingPoint[3]=0.898;
  }
  else if(Tagger=="combinedSecondaryVertexV1BJetTags") {
    WorkingPoint[1]=0.405;
    WorkingPoint[2]=0.783;
    WorkingPoint[3]=0.920;
  }
  else if(Tagger=="combinedSecondaryVertexSoftPFLeptonV1BJetTags") {
    WorkingPoint[1]=0.527;
    WorkingPoint[2]=0.756;
    WorkingPoint[3]=0.859;
  }
  WorkingPoint[4]=1.;
  
  for(int f=0; f<3; f++) for(int i=0; i<5; i++) Integral[f][i]=Shape[f]->Integral(Shape[f]->FindBin(WorkingPoint[i]), Shape[f]->GetNbinsX()+1);
  
  std::cout << " - BTagInterface initialized, tagger: " << bTag << std::endl;
}

BTagInterface::~BTagInterface() {
  BDFile->Close();
  ShapesFile->Close();
}





void BTagInterface::FillBTagVector(const edm::Event& iEvent, std::vector<pat::Jet>& Vect) {
  bool isMC(!iEvent.isRealData());
  // Loop on jet vector
  for(unsigned int i=0; i<Vect.size(); i++) {
    float Reshaped(0.), ReshapedUp(0.), ReshapedDown(0.);
    if(!isMC) Reshaped=ReshapedUp=ReshapedDown=Vect[i].bDiscriminator(Tagger);
    else { // Reshape discriminator only for MC
      Reshaped=GetReshapedDiscr(Vect[i], 0);
      ReshapedUp=GetReshapedDiscr(Vect[i], +1);
      ReshapedDown=GetReshapedDiscr(Vect[i], -1);
    }
    // Write discriminators
    Vect[i].addBDiscriminatorPair(std::pair<std::string, float>(std::string(Tagger+"Reshaped"), Reshaped));
    Vect[i].addBDiscriminatorPair(std::pair<std::string, float>(std::string(Tagger+"ReshapedUp"), ReshapedUp));
    Vect[i].addBDiscriminatorPair(std::pair<std::string, float>(std::string(Tagger+"ReshapedDown"), ReshapedDown));
  }
}

float BTagInterface::GetScaleFactor(int wp, float x) {
  if(wp<=0 || wp>=4) return 1.;
  if(Tagger=="combinedSecondaryVertexBJetTags") {
    if(wp==1) return 0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)));
    if(wp==2) return (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
    if(wp==3) return (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
  }
  else if(Tagger=="combinedSecondaryVertexV1BJetTags") {
    if(wp==1) return 1.7586*((1.+(0.799078*x))/(1.+(1.44245*x)));
    if(wp==2) return 0.952067+(-2.00037e-05*x);
    if(wp==3) return (0.912578+(0.000115164*x))+(-2.24429e-07*(x*x));
  }
  else if(Tagger=="combinedSecondaryVertexSoftPFLeptonV1BJetTags") {
    if(wp==1) return 0.970168*((1.+(0.00266812*x))/(1.+(0.00250852*x)));
    if(wp==2) return ((0.939238+(0.000278928*x))+(-7.49693e-07*(x*x)))+(2.04822e-10*(x*(x*x)));
    if(wp==3) return (0.928257+(9.3526e-05*x))+(-4.1568e-07*(x*x));
  }
  return 1.;
}

float BTagInterface::GetScaleFactor(int wp, float pt, float eta, int sigma) {
  if(!isBDFile) return 1.;
  eta=fabs(eta);
  if(wp<=0 || wp>=4) return 1.;
  
  //float sf=SFb[wp-1]->Interpolate(pt, eta) + sigma*SFb[wp-1]->GetBinError(SFb[wp-1]->FindBin(pt, eta));
  float sf=SFb[wp-1]->GetBinContent(SFb[wp-1]->FindBin(pt, eta)) + sigma*SFb[wp-1]->GetBinError(SFb[wp-1]->FindBin(pt, eta));
  if(sf!=sf || isinf(sf)) return -1.;
  return sf;
}

float BTagInterface::GetScaleFactorMistag(int wp, float pt, float eta, int sigma) {
  if(!isBDFile) return 1.;
  eta=fabs(eta);
  if(wp<=0 || wp>=4) return 1.;
  
  //float sf=SFl[wp-1]->Interpolate(pt, eta) + sigma*SFl[wp-1]->GetBinError(SFl[wp-1]->FindBin(pt, eta));
  float sf=SFl[wp-1]->GetBinContent(SFl[wp-1]->FindBin(pt, eta)) + sigma*SFl[wp-1]->GetBinError(SFl[wp-1]->FindBin(pt, eta));
  if(sf!=sf || isinf(sf)) return -1.;
  return sf;
}


float BTagInterface::GetWorkingPoint(int wp) {
  if(wp<0 || wp>4) {
    std::cout << " - BTagInterface Error: working point not defined" << std::endl;
    return -1.;
  }
  return WorkingPoint[wp];
}

float BTagInterface::GetNewWorkingPoint(int fl, int wp, float pt, float eta, int sigma) {
  if(fl<0 || fl>2) {
    std::cout << " - BTagInterface Error: unrecognized flavour" << std::endl;
    return -1.;
  }
  if(wp<0 || wp>4) {
    std::cout << " - BTagInterface Error: working point not defined" << std::endl;
    return -1.;
  }
  if(wp==0 || wp==4) return WorkingPoint[wp];
  if(!isShapesFile) return WorkingPoint[wp];
  // Get original integral
  //float Int=Shape[fl]->Integral(Shape[fl]->FindBin(WorkingPoint[wp]), Shape[fl]->GetNbinsX()+1);
  float Int=Integral[fl][wp]; // Saves time
  // Scale integral by Scale Factor with Error
  if(fl==0) Int/=GetScaleFactor(wp, pt, eta, sigma);
  else if(fl==1) Int/=GetScaleFactor(wp, pt, eta, 2*sigma);
  else Int/=GetScaleFactorMistag(wp, pt, eta, sigma);
  // Find new Discriminator value
  // Fast scan, start from wp and go up
  for(int i=Shape[fl]->GetNbinsX()+1; i>0; i=i-100) {
    if(Shape[fl]->Integral(i-100, Shape[fl]->GetNbinsX()+1)>=Int) {
      for(; i>0; i=i-10) {
        if(Shape[fl]->Integral(i-10, Shape[fl]->GetNbinsX()+1)>=Int) {
          for(; i>0; i=i-1) {
            if(Shape[fl]->Integral(i, Shape[fl]->GetNbinsX()+1)>=Int) {
              return ((float)i-0.5)/Shape[fl]->GetNbinsX();
            }
          }
        }
      }
    }
  }
  // Integral Interpolation
//  for(int i=Shape[fl]->GetNbinsX()+1; i>0; i--) {
//    if(Shape[fl]->Integral(i, Shape[fl]->GetNbinsX()+1)>=Int) {
//      float x0=Shape[fl]->Integral(i, Shape[fl]->GetNbinsX()+1);
//      float x1=Shape[fl]->Integral(i+1, Shape[fl]->GetNbinsX()+1);
//      float y0=Shape[fl]->GetXaxis()->GetBinLowEdge(i);
//      float y1=Shape[fl]->GetXaxis()->GetBinLowEdge(i+1);
//      return y0 + (Int-x0)*((y1-y0)/(x1-x0));
//    }
//  }
  //std::cout << " - BTagWeight Warning: new working point " << wp << " not found" << std::endl;
  return WorkingPoint[wp];
}

float BTagInterface::GetReshapedDiscr(pat::Jet& jet, int sigma) {
  // Flavour
  int fl(-1);
  if(abs(jet.partonFlavour())==5) fl=0;
  if(abs(jet.partonFlavour())==4) fl=1;
  if(abs(jet.partonFlavour())<4 || abs(jet.partonFlavour())==21) fl=2;
  // Discrminator
  float discr=jet.bDiscriminator(Tagger);
  if(discr<0.01 || discr>0.99) return discr;
  if(fl<0 || fl>2) return discr;
  
  // Find boundary old Discr values
  int i0(0), i1(4);
  float x0(0.), x1(1.);
  for(int i=1; i<5; i++) {
    if(discr<=WorkingPoint[i]) {
      x0=WorkingPoint[i-1];
      x1=WorkingPoint[i];
      i0=i-1;
      i1=i;
      break;
    }
  }
  // Find boundary new Discr values
  float y0=GetNewWorkingPoint(fl, i0, jet.pt(), jet.eta(), sigma);
  float y1=GetNewWorkingPoint(fl, i1, jet.pt(), jet.eta(), sigma);
  // Interpolate
  return y0 + (discr-x0)*((y1-y0)/(x1-x0));
}

