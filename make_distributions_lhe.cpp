// c++ -o checkMomentum_00 `root-config --glibs --cflags` -lm checkMomentum_00.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "TH1.h"
#include "TFile.h"

#include "TLorentzVector.h"

// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;
using namespace std ;

TLorentzVector buildP (const LHEF::HEPEUP & event, int iPart)
{
  TLorentzVector dummy ;
  dummy.SetPxPyPzE (
      event.PUP.at (iPart).at (0), // px
      event.PUP.at (iPart).at (1), // py
      event.PUP.at (iPart).at (2), // pz
      event.PUP.at (iPart).at (3) // E
    ) ;
  return dummy ;  
}


int main(int argc, char ** argv) 
{
  if(argc < 2)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe " << endl ;
      return -1;
    }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  TH1F dielectron_mass ("dielectron_mass", "dielectron_mass", 100, 0, 300) ;
  TH1F diquark_mass ("diquark_mass", "diquark_mass", 100, 0, 100) ;
  TH1F leading_lep_pt ("leading_lep_pt", "leading_lep_pt", 100, 0, 200) ;
  TH1F trailing_lep_pt ("trailing_lep_pt", "trailing_lep_pt", 100, 0, 200) ;
  TH1F leading_quark_pt ("leading_quark_pt", "leading_quark_pt", 100, 0, 400) ;
  TH1F trailing_quark_pt ("trailing_quark_pt", "trailing_quark_pt", 100, 0, 200) ;

  std::vector<TLorentzVector> electrons;

  //PG loop over input events
  while (reader.readEvent ()) 
    {

      std::vector<TLorentzVector> electrons;
      std::vector<TLorentzVector> quarks;

      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
           // outgoing particles          
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
               if (abs (reader.hepeup.IDUP.at (iPart)) == 11) 
                 {     
                   TLorentzVector vec = buildP (reader.hepeup, iPart) ;
		   electrons.push_back(vec);
                 }

	       if(abs (reader.hepeup.IDUP.at (iPart)) == 5 || abs (reader.hepeup.IDUP.at (iPart)) == 6)
		 std::cout << "abs (reader.hepeup.IDUP.at (iPart)) = " << abs (reader.hepeup.IDUP.at (iPart)) << std::endl;

               if (abs (reader.hepeup.IDUP.at (iPart)) == 1 || abs (reader.hepeup.IDUP.at (iPart)) == 2 || abs (reader.hepeup.IDUP.at (iPart)) == 3 || abs (reader.hepeup.IDUP.at (iPart)) == 4)
                 {
                   TLorentzVector vec = buildP (reader.hepeup, iPart) ;
                   quarks.push_back(vec);
                 }
 
             } // outgoing particles
        } // loop over particles in the event

      assert(quarks.size() == 2);
      assert(electrons.size() == 2);

      dielectron_mass.Fill( (electrons[0] + electrons[1]).M() );
      diquark_mass.Fill( (quarks[0] + quarks[1]).M() );
      if (electrons[0].Pt() > electrons[1].Pt()){
	leading_lep_pt.Fill(electrons[0].Pt());
	trailing_lep_pt.Fill(electrons[1].Pt());
      }
      else {
        leading_lep_pt.Fill(electrons[1].Pt());
        trailing_lep_pt.Fill(electrons[0].Pt());
      }

      if (quarks[0].Pt() > quarks[1].Pt()){
        leading_quark_pt.Fill(quarks[0].Pt());
        trailing_quark_pt.Fill(quarks[1].Pt());
      }
      else {
        leading_quark_pt.Fill(quarks[1].Pt());
        trailing_quark_pt.Fill(quarks[0].Pt());
      }

    } //PG loop over input events

  TFile f ("output_distributions_00.root", "recreate") ;
  dielectron_mass.Write();
  diquark_mass.Write();
  leading_lep_pt.Write();
  trailing_lep_pt.Write();
  leading_quark_pt.Write();
  trailing_quark_pt.Write();
  f.Close () ;

  return 0 ;
}
