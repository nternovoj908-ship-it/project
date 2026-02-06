// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file TimepixDetector.cc
/// \brief Implementation of the B1::TimepixDetector class

#include "RootManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4VSensitiveDetector.hh"
#include "TimepixDetector.hh"
#include "TimepixHit.hh"
#include "G4THitsCollection.hh"
#include "G4VHitsCollection.hh"

#include "G4SDManager.hh"

using namespace std;

TimepixDetector::TimepixDetector(G4String name) 
: G4VSensitiveDetector(name)
{ 
	collectionName.insert("TimepixHitCollection");
}

TimepixDetector::~TimepixDetector()
{ }

G4bool TimepixDetector::ProcessHits(G4Step* step, G4TouchableHistory* hist)
{ 
	// Создаём hit через коллекцию — Geant4 сама управляет памятью
	auto* HC = THC;
	if (!HC) return false;

	TimepixHit* hit = new TimepixHit();

	G4StepPoint* point = step->GetPreStepPoint();
	G4ThreeVector pos = point->GetPosition();

	G4double edep = step -> GetTotalEnergyDeposit();
	G4double kin = step -> GetTrack() -> GetKineticEnergy();
	G4double fullkin = step -> GetTrack() -> GetKineticEnergy();
	
	hit -> AddX(pos.x());
	hit -> AddY(pos.y());
	hit -> AddE(edep);
	hit -> AddK(kin);
	hit -> AddFE(fullkin);

	// Добавляем в коллекцию — Geant4 сама удалит
	HC -> insert(hit);

	return true;
}

void TimepixDetector::Initialize(G4HCofThisEvent* HCE)
{ 
	THC = new TimepixHitCollection(GetName(), collectionName[0]);
	static G4int HCID = -1;
	if (HCID < 0 ) HCID = GetCollectionID(0);
	HCE -> AddHitsCollection(HCID, THC);
}

void TimepixDetector::EndOfEvent(G4HCofThisEvent* HCE)
{ 
	G4int id = G4SDManager::GetSDMpointer() -> GetCollectionID("Timepix1/TimepixHitCollection");
	auto* coll = (TimepixHitCollection*)HCE -> GetHC(id);	

	G4int number = coll -> entries();

	if (number == 0) return;

	G4double hit_x = 0;
	G4double hit_y = 0;
	G4double energy = 0;
	G4double kinetic = 0;
	G4double fullkin = 0;

	for (int i = 0; i < number; i++)
	{
		auto* hit = (*coll)[i];

		hit_x+= hit -> GetX();
		hit_y+= hit -> GetY();
		energy+= hit -> GetE();
		kinetic+= hit -> GetK();	
		fullkin+= hit -> GetFE();	
	}

	energy = energy/keV;
	hit_x = hit_x/mm/number;
	hit_y = hit_y/mm/number;
	kinetic = kinetic/keV/number;
	fullkin = fullkin/keV/number;
	
	MyROOTManager* RM = MyROOTManager::GetPointer();		
	RM -> FillHist(hit_x, hit_y);
	RM -> FillTree(hit_x, hit_y, energy, kinetic, fullkin);

	// ---> Ручная очистка коллекции после обработки
	THC->clear();  // <-- Это гарантирует, что все TimepixHit* будут удалены
}
