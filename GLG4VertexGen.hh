////////////////////////////////////////////////////////////////////
// Last svn revision: $Id$ 
////////////////////////////////////////////////////////////////////
#ifndef __GLG4VertexGen_h__
#define __GLG4VertexGen_h__ 1
/** @file
 Declares GenericLAND global vertex generator classes for primary events,
 (See note on GenericLAND generators for more information.)

 This file is part of the GenericLAND software library.
 $Id$

 @author G.Horton-Smith, August 3, 2001
 Aleksandra Bialek: possiblity to vary number of particles per event:
 additional parameter fRandom.
*/

#include "globals.hh"
#include <stdio.h> // for FILE
#include <G4ThreeVector.hh>

class G4Event;
class G4Track;
class G4PrimaryVertex;
class G4ParticleDefinition;
class GLG4PrimaryGeneratorAction;

/** Virtual base class for vertex generators */
class GLG4VertexGen {
public:
  GLG4VertexGen(const char *arg_dbname="vertex") : fDbname(arg_dbname) { };
  virtual ~GLG4VertexGen() { };
  // Generate and add new vertex to this event.  Position and
  // time of vertex are offset from 0 by dx and dt.  (These
  // are usually derived from GLG4PosGen and GLG4TimeGen.)
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
				      G4ThreeVector &dx,
				      G4double dt)=0;
  // adds vertices.
  virtual void SetState( G4String newValues ) = 0;
  // sets filename or other information needed by vertex generator
  virtual G4String GetState() = 0;
  // returns the current state information in a form that can be understood
  // by SetState (and, hopefully, a well-informed human)
  /** This method tells us that you can't limit the output energies (overload if true)*/
  virtual bool ELimitable() { return false; };
  /** Optional method to limit the energies, must be supplied if ELimitable returns true */
  virtual void LimitEnergies( float, float ) { };
  /** Optional method to return maximum energy possible for generator */
  virtual float EMaximum() { return 99999; };
  /** Optional method to return minimum energy possible for generator */
  virtual float EMinimum() { return 0; };
  
protected:
  G4String fDbname; // used for GLG4param key prefix
};


/** vertex generator that can generate a primary vertex with one more
    particles of a given type, direction, energy, and polarization.
    Allows for randomly isotropic direction and random transverse polarization
    of spin-1, mass=0 particles */
class GLG4VertexGen_Gun : public GLG4VertexGen {
public:
  GLG4VertexGen_Gun(const char *arg_dbname="gun");
  virtual ~GLG4VertexGen_Gun();
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
				      G4ThreeVector &dx,
				      G4double dt);
  // generates a primary vertex with given particle type, direction, energy,
  // and consistent polarization.

  virtual void SetState( G4String newValues );
  // format: particle_name  dir_x dir_y dir_z  kinetic_energy  polx poly polz
  // If dir_x==dir_y==dir_z==0, the directions are isotropic.
  // If particle has mass==0 and spin==1, final polarization will be
  // projected into plane perpindicular to momentum and made a unit vector;
  // if polarization has zero magnitude, a polarization is chosen randomly.

  virtual G4String GetState();
  // returns current state formatted as above

public:
  // the following useful static const data should be universally accessable
  // (I copied it from the G4IonTable source code, where it is privatized
  // with no accessor functions.)
  enum { fNumberOfElements = 110};
  static const char * fTheElementNames[fNumberOfElements];

protected:
  G4ParticleDefinition *fPDef;
  G4ThreeVector fMom;
  G4double fKe;
  G4ThreeVector fPol;
  G4int fMultiplicity;//! small modification to have random number of initial particles: 
  G4int fRandom; //! fRandom==1 than:max = fMultiplicity*G4UniformRand()+1;
  //! default: fixed number of particles (Aleksandra)
};

/////////////////////////////   GUN2
/** vertex generator that can generate a primary vertex with one more
    particles of a given type, direction, energy, and polarization.
    Allows for randomly isotropic direction and random transverse polarization
    of spin-1, mass=0 particles */
class GLG4VertexGen_Gun2 : public GLG4VertexGen {
public:
  GLG4VertexGen_Gun2(const char *arg_dbname="gun2");
  virtual ~GLG4VertexGen_Gun2();
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
				      G4ThreeVector &dx,
				      G4double dt);
  // generates a primary vertex with given particle type, direction, energy,
  // and consistent polarization.

  virtual void SetState( G4String newValues );
  // format: particle_name  dir_x dir_y dir_z  kinetic_energy  polx poly polz
  // If dir_x==dir_y==dir_z==0, the directions are isotropic.
  // If particle has mass==0 and spin==1, final polarization will be
  // projected into plane perpindicular to momentum and made a unit vector;
  // if polarization has zero magnitude, a polarization is chosen randomly.

  virtual G4String GetState();
  // returns current state formatted as above

public:
  // the following useful static const data should be universally accessable
  // (I copied it from the G4IonTable source code, where it is privatized
  // with no accessor functions.)
  enum { fNumberOfElements = 110};
  static const char * fTheElementNames[fNumberOfElements];

protected:
  G4ParticleDefinition *fPDef;
  G4ThreeVector fMom;
  G4double fMomTheta;
  G4double fKe1;
  G4double fKe2;
  G4ThreeVector fPol;
  G4int fMultiplicity;
  G4int fRandom; 
};

/////////////////////////////   GUN3
/** vertex generator that can generate a primary vertex with one more
    particles of a given type, direction, energy, and polarization, and time.
    Allows for randomly isotropic direction and random transverse polarization
    of spin-1, mass=0 particles */
class GLG4VertexGen_Gun3 : public GLG4VertexGen {
public:
  GLG4VertexGen_Gun3(const char *arg_dbname="gun3");
  virtual ~GLG4VertexGen_Gun3();
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
				      G4ThreeVector &dx,
				      G4double dt);
  // generates a primary vertex with given particle type, direction, energy,
  // and consistent polarization.

  virtual void SetState( G4String newValues );
  // format: particle_name  dir_x dir_y dir_z  kinetic_energy  polx poly polz
  // If dir_x==dir_y==dir_z==0, the directions are isotropic.
  // If particle has mass==0 and spin==1, final polarization will be
  // projected into plane perpindicular to momentum and made a unit vector;
  // if polarization has zero magnitude, a polarization is chosen randomly.

  virtual G4String GetState();
  // returns current state formatted as above

public:
  // the following useful static const data should be universally accessable
  // (I copied it from the G4IonTable source code, where it is privatized
  // with no accessor functions.)
  enum { fNumberOfElements = 110};
  static const char * fTheElementNames[fNumberOfElements];

protected:
  G4ParticleDefinition *fPDef;
  G4ThreeVector fMom;
  G4double fMomTheta;
  G4double fKe1;
  G4double fKe2;
  G4ThreeVector fPol;
  G4int fMultiplicity;
  G4int fRandom; 
  G4double fTime;
};


/** vertex generator that generates a primary vertex based on
    information from an ASCII stream.
    The ASCII stream contains information in a HEPEVT-like style.
*/
class GLG4VertexGen_HEPEvt : public GLG4VertexGen {
public:
  GLG4VertexGen_HEPEvt(const char *arg_dbname);
  virtual ~GLG4VertexGen_HEPEvt();
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
				      G4ThreeVector &dx,
				      G4double dt)=0;

  // Generates a primary vertex based on information from an ASCII stream.
  // The ASCII stream contains information in a HEPEVT style.

  virtual void SetState( G4String newValues );
  // The argument is a Perl "open" style filename argument:
  // if the filename ends with '|', the filename is interpreted as
  // a command which pipes output to us; otherwise, if the filename begins
  // with '<' or an ordinary character, the file is opened for input.
  // The "pipe" style can be used to read a gzip'ped file or to start
  // a program which feeds events to us directly.
  virtual G4String GetState();
  // returns current state formatted as above

  void Open(const char *argFilename);
  void GetDataLine(char *buffer, size_t size);
  void Close();

  enum { kIonCodeOffset = 9800000,       // nuclei have codes like 98zzaaa
	 kPDGcodeModulus=10000000,       // PDG codes are 7 digits long
     kISTHEP_ParticleForTracking=1,  // only ISTHEP==1 are for tracking
	 kISTHEP_InformatonMin=100,      // 100 and above are "informatons"
	 kISTHEP_Max=213                 // <= MAXINT / kPDGcodeModulus
  };
  
protected:
  G4String fFileName;
  FILE *   fFile;
  bool fIsPipe;
};



#endif
