#include "TimepixHit.hh"
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"

// Определяем allocator (важно!)
G4ThreadLocal G4Allocator<TimepixHit>* TimepixHitAllocator = 0;

TimepixHit::TimepixHit() : G4VHit()
{ }

TimepixHit::~TimepixHit()
{ }

void TimepixHit::Print()
{ }

void TimepixHit::Draw()
{ }

void TimepixHit::AddX(G4double a)
{
	pos_x = a;	
}

void TimepixHit::AddY(G4double a)
{
	pos_y = a;	
}

void TimepixHit::AddE(G4double a)
{
	edep = a;	
}

void TimepixHit::AddK(G4double a)
{
	ekin = a;	
}

void TimepixHit::AddFE(G4double a)
{
	fullenergy = a;	
}

G4double TimepixHit::GetX()
{
	return pos_x;
}

G4double TimepixHit::GetY()
{
	return pos_y;
}

G4double TimepixHit::GetE()
{
	return edep;
}

G4double TimepixHit::GetK()
{
	return ekin;
}

G4double TimepixHit::GetFE()
{
	return fullenergy;
}

