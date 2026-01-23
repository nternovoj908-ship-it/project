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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh" // <-- Добавлено для цилиндрических форм
#include <cstdio>    // <-- Добавлено для работы с файлом

// Добавляем подключения для чувствительного детектора
#include "TimepixDetector.hh"
#include "G4SDManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20 * cm, env_sizeZ = 60 * cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2 * env_sizeXY;
  G4double world_sizeZ = 1.2 * env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",  // its name
                            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                      env_mat,  // its material
                                      "Envelope");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicEnv,  // its logical volume
                    "Envelope",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking

  //
  // Collimator (replacing Shape1)
  //

  G4Material* collimator_mat = nist->FindOrBuildMaterial("G4_W"); // Вольфрам
  G4ThreeVector collimator_pos = G4ThreeVector(0, 0, -7 * cm); // Та же Z-координина, что была у Shape1

  G4double collimator_outer_radius = 2.203869 * cm;
  G4double collimator_inner_radius = 0 * cm;
  G4double collimator_thickness = 1. * mm;
  G4double collimator_phi1 = 0 * deg;
  G4double collimator_phi2 = 360.0 * deg;

  auto solidCollimator = new G4Tubs("Collimator",
                                    collimator_inner_radius,
                                    collimator_outer_radius,
                                    collimator_thickness / 2.,
                                    collimator_phi1,
                                    collimator_phi2);

  auto logicCollimator = new G4LogicalVolume(solidCollimator,
                                             collimator_mat,
                                             "LogicCollimator");

  auto physCollimator = new G4PVPlacement(nullptr, // no rotation
                                          collimator_pos,
                                          logicCollimator,
                                          "PhysCollimator",
                                          logicEnv, // размещаем внутри "Envelope"
                                          false,   // no boolean operation
                                          0,       // copy number
                                          checkOverlaps);

  //
  // Hole placement in collimator (based on mask file)
  //

  FILE* mask_conf = fopen("../mask-mura-61.txt", "r");
  if (!mask_conf) {
      G4Exception("DetectorConstruction::Construct()", "NoMaskFile", FatalException,
                  "Could not open mask file ../mask-mura-61.txt");
  }

  const int mask_size = 61;
  int** mask_pattern = new int*[mask_size];

  for (int i = 0; i < mask_size; i++) {
      mask_pattern[i] = new int[mask_size];
      for (int j = 0; j < mask_size; j++) {
          int tmp;
          fscanf(mask_conf, "%d", &tmp);
          mask_pattern[i][j] = tmp;
      }
  }

  fclose(mask_conf);

  //
  // Hole parameters
  //

  G4Material* hole_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4double hole_outer_radius = 0.015 * cm;
  G4double hole_thickness = 1. * mm;
  G4double hole_step = 0.036129 * cm;

  auto solidHole = new G4Tubs("Hole",
                              0, // inner radius
                              hole_outer_radius,
                              hole_thickness / 2.,
                              0 * deg,
                              360. * deg);

  auto logicHole = new G4LogicalVolume(solidHole,
                                       hole_mat,
                                       "LogicHole");

  //
  // Place holes according to mask pattern
  //

  int Inversion = 0; // 0 - off, 1 - on (определяет, где дыры)

  for (int i = 0; i < mask_size; i++) {
      for (int j = 0; j < mask_size; j++) {
          G4double hole_pos_x = (j - (mask_size - 1) / 2.0) * hole_step;
          G4double hole_pos_y = ((mask_size - 1) / 2.0 - i) * hole_step;

          G4ThreeVector hole_pos = G4ThreeVector(hole_pos_x, hole_pos_y, 0);

          if (mask_pattern[i][j] == Inversion) {
              new G4PVPlacement(nullptr, // no rotation
                                hole_pos,
                                logicHole,
                                "PhysHole",
                                logicCollimator, // размещаем в логическом объёме коллиматора
                                false,
                                i * mask_size + j,
                                checkOverlaps);
          }
      }
  }

  //
  // Free allocated memory
  //

  for (int i = 0; i < mask_size; i++) {
      delete[] mask_pattern[i];
  }
  delete[] mask_pattern;

  //
  // New detector (from second file)
  //

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
  G4ThreeVector detector_pos = G4ThreeVector(0, 0, 4.4216 * cm); // Позиция нового детектора

  G4double detector_width = 14.08 * mm;
  G4double detector_thickness = 1 * mm;

  auto solidDetector = new G4Box("Detector",
                                 detector_width / 2.,
                                 detector_width / 2.,
                                 detector_thickness / 2.);

  auto logicDetector = new G4LogicalVolume(solidDetector,
                                           detector_mat,
                                           "LogicDetector");

  new G4PVPlacement(nullptr, // no rotation
                    detector_pos,
                    logicDetector,
                    "PhysDetector",
                    logicEnv, // размещаем внутри "Envelope"
                    false,
                    0,
                    checkOverlaps);

  // Set detector as scoring volume
  fScoringVolume = logicDetector;

  // Set up sensitive detector
  TimepixDetector* timepix = new TimepixDetector("/Timepix1");
  G4SDManager* SDM = G4SDManager::GetSDMpointer();
  SDM->AddNewDetector(timepix);
  logicDetector->SetSensitiveDetector(timepix);

  //
  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
