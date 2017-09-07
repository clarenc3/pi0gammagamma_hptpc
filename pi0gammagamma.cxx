#include "pi0gammagamma.h"

int main(int argc, char* argv[]) {

  if (argc != 2) {
    Usage();
    exit(-1);
  }

  pi0gammagamma(std::string(argv[1]));

  return 0;
};

void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "./pi0gammagamma.exe ROOT_FILE" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC0pi, 1 = CC1pi, 2 = CC1pi0" << std::endl;
  std::cout << "MAX_MOM is the rough momentum scale of the neutrinos" << std::endl;

  return;
}

void pi0gammagamma(std::string fileName) {

  TFile *f = TFile::Open((fileName).c_str(),"open");
  TTree *tn = (TTree*)(f->Get("neuttree"));
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

  long int nevents = tn->GetEntries();

  int eventCnt = 0;

  double thetaHigh = -999;
  double pgammaHigh = -999;
  double thetaLow = -999;
  double pgammaLow = -999;
  double thetagamma = -999;
  int nPi0 = -999;
  double pi0mom = -999;
  double pi0th = -999;

  TFile *output = new TFile((fileName+"_output.root").c_str(), "recreate");
  TTree *tree = new TTree("tree", "tree");
  tree->Branch("thgammaHigh", &thetaHigh);
  tree->Branch("pgammaHigh", &pgammaHigh);
  tree->Branch("thgammaLow", &thetaLow);
  tree->Branch("pgammaLow", &pgammaLow);
  tree->Branch("thgammaHighgammaLow", &thetagamma);
  tree->Branch("ppi0", &pi0mom);
  tree->Branch("thpi0", &pi0th);
  tree->Branch("nPi0", &nPi0);

  std::cout << "Number of events: " << nevents << std::endl;

  // 20 messages will do the job
  int countwidth = double(nevents)/double(20);

  // Random number for generator the theta and phi in pi0 rest-frame
  TRandom3 *rand = new TRandom3();

  for (int j = 0; j < nevents; j++) {

    // Reset the variables
    thetaHigh = -999;
    pgammaHigh = -999;
    thetaLow = -999;
    pgammaLow = -999;
    thetagamma = -999;
    nPi0 = 0;
    pi0mom = -999;
    pi0th = -999;

    tn->GetEntry(j);

    if (j % countwidth == 0) {
      std::cout << "On event " << j << "/" << nevents << " (" << int(double(j)/double(nevents)*100.0) << "%)" << std::endl;
    }

    TLorentzVector Ppi0;

    // Loop over the particle stack, find the highest momentum pi0
    for (int k = 2; k < nvect->Npart(); ++k) {
      int PID = (nvect->PartInfo(k))->fPID;
      if (PID == 111) {
        nPi0++;
        // Pick the highest momentum particle
        if (nvect->PartInfo(k)->fP.Mag() > Ppi0.Mag()) {
          Ppi0 = nvect->PartInfo(k)->fP;
        }
      }
    }
    // Skip events without pi0
    if (nPi0 == 0) continue;
    // Increase good event counter
    eventCnt++;

    // Make the angles relative incoming neutrino beam axis
    TLorentzVector Pnu = nvect->PartInfo(0)->fP;

    // Calculate pi0 momentum
    pi0mom = Ppi0.Vect().Mag();
    pi0th = 180.0/M_PI*Ppi0.Vect().Angle(Pnu.Vect());

    // Make the gammas at random angles in the pi0 rest-frame with momentum balanced
    double mass = Ppi0.Mag();
    // The energy of each gamma
    double energy = mass/2.0;
    // Generate isotropically: our only constraint is the energy of the system
    double theta = M_PI*rand->Rndm();
    double phi = 2*M_PI*rand->Rndm();

    // Make the Cartesian x,y,z
    double x = energy*sin(theta)*cos(phi);
    double y = energy*sin(theta)*sin(phi);
    double z = energy*cos(theta);

    // Now make the photons
    TLorentzVector gamma1(x, y, z, energy);
    TLorentzVector gamma2(-x, -y, -z, energy);

    // Boost out of pi0 rest-frame into lab frame
    gamma1.Boost(Ppi0.BoostVector());
    gamma2.Boost(Ppi0.BoostVector());

    // Select highest momentum gamma to be "high" and lowest momentum to be "low"
    if (gamma1.Vect().Mag() > gamma2.Vect().Mag()) {
      thetaHigh = 180.0/M_PI*gamma1.Vect().Angle(Pnu.Vect());
      pgammaHigh = gamma1.Vect().Mag();
      thetagamma = 180.0/M_PI*gamma1.Vect().Angle(gamma2.Vect());
      thetaLow = 180.0/M_PI*gamma2.Vect().Angle(Pnu.Vect());
      pgammaLow = gamma2.Vect().Mag();
    } else {
      thetaHigh = 180.0/M_PI*gamma2.Vect().Angle(Pnu.Vect());
      pgammaHigh = gamma2.Vect().Mag();
      thetagamma = 180.0/M_PI*gamma2.Vect().Angle(gamma1.Vect());
      thetaLow = 180.0/M_PI*gamma1.Vect().Angle(Pnu.Vect());
      pgammaLow = gamma1.Vect().Mag();
    }

    tree->Fill();

  } // End for loop over all events

  std::cout << eventCnt << "/" << nevents << " (" << int(double(eventCnt)/double(nevents)*100.0) << "%) events passed the cuts" << std::endl;
  output->cd();
  tree->Write();
  output->Close();
  std::cout << "Wrote to " << output->GetName() << std::endl;

};
