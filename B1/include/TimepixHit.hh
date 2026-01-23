#ifndef TimepixHit_h
#define TimepixHit_h 1

#include "G4THitsCollection.hh"
#include "globals.hh"
#include "G4VHit.hh"


//A class inherited from default to make G4SensitiveDetector work

class TimepixHit : public G4VHit
{
	public:
	TimepixHit();
	~TimepixHit();

	void Print();
	void Draw();

	void AddX (G4double);
	void AddY (G4double);
	void AddE (G4double);
	void AddK (G4double);
	void AddFE (G4double);
	G4double GetX();
	G4double GetY();
	G4double GetE();
	G4double GetK();
	G4double GetFE();


	private:
	G4double pos_x, pos_y, edep, ekin, fullenergy;	

};



#endif

typedef G4THitsCollection <TimepixHit> TimepixHitCollection;
