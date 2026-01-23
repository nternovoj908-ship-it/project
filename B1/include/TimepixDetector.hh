#ifndef TimepixDetector_h
#define TimepixDetector_h 1

#include "G4THitsCollection.hh"
#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "TimepixHit.hh"

//A class inherited from default to make G4SensitiveDetector work

class TimepixDetector : public G4VSensitiveDetector
{
	public:
	TimepixDetector(G4String);
	~TimepixDetector();
	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	void Initialize(G4HCofThisEvent*);
	void EndOfEvent(G4HCofThisEvent*);

	private:
	TimepixHitCollection* THC;

};

typedef G4THitsCollection <TimepixHit> TimepixHitCollection;
#endif
