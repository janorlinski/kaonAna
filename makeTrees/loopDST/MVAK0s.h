#include "TMVA/Reader.h"

TString MVAWeightfileK0s[4] = {"/lustre/nyx/hades/user/lchlad/flow/scripts/myTMVA/training/dataset/weights/TMVAClassification_400000_events_Jan20_MPC2_WKM_MLP.weights.xml",
                               "/lustre/nyx/hades/user/lchlad/flow/scripts/myTMVA/training/dataset/weights/TMVAClassification_400000_events_Jan20_SPC_WKM_MLP.weights.xml",
                               "/lustre/nyx/hades/user/lchlad/flow/scripts/myTMVA/training/dataset/weights/TMVAClassification_200000_events_MPC1_WKM_MLP.weights.xml",
                               "/lustre/nyx/hades/user/lchlad/flow/scripts/myTMVA/training/dataset/weights/TMVAClassification_400000_events_Jan20_MPC2_WKM_MLP.weights.xml"
                              };


TMVA::Reader* MVAReaderK0s;

void InitMVAK0s(Float_t* TMVAd1, Float_t* TMVAd2, Float_t* TMVAd3, Float_t* TMVAdVer, Float_t* TMVAdMin, Float_t* kaonMom, Int_t iPrecuts) {
    MVAReaderK0s = new TMVA::Reader("!V:Color:!Silent");

    MVAReaderK0s->AddVariable("d1",   TMVAd1);
    MVAReaderK0s->AddVariable("d2",   TMVAd2);
    MVAReaderK0s->AddVariable("d3",   TMVAd3);
    MVAReaderK0s->AddVariable("dVer", TMVAdVer);
    MVAReaderK0s->AddVariable("dMin", TMVAdMin);
    if( MVAWeightfileK0s[iPrecuts].Contains("WKM") ) MVAReaderK0s->AddVariable("kaonMom",  kaonMom);

    MVAReaderK0s->BookMVA("MLPK0s", MVAWeightfileK0s[iPrecuts]);
}

Double_t EvalMVAK0s() { return MVAReaderK0s->EvaluateMVA("MLPK0s"); }
