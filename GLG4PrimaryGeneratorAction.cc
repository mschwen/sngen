////////////////////////////////////////////////////////////////////
// Last svn revision: $Id$ 
////////////////////////////////////////////////////////////////////
// This file is part of the GenericLAND software library.
// $Id$
//
// new GLG4PrimaryGeneratorAction.cc by Glenn Horton-Smith, August 3-17, 2001

////////////////////////////////////////////////////////////////
// GLG4PrimaryGeneratorAction
////////////////////////////////////////////////////////////////

#include <RAT/GLG4PrimaryGeneratorAction.hh>
#include <RAT/GLG4PrimaryGeneratorMessenger.hh>

#include "globals.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include <G4RunManager.hh>
#include "CLHEP/Units/SystemOfUnits.h"
#include "Randomize.hh"

#include <RAT/GLG4Gen.hh>
#include <RAT/GLG4TimeGen.hh>
#include <RAT/GLG4VertexGen.hh>            // for vertex generator
#include <RAT/GLG4param.hh>                 // for GLG4param
#include <RAT/Factory.hh>
#include <RAT/GLG4PosGen.hh>
#include <RAT/PosGen_Point.hh>
#include <RAT/PosGen_Fill.hh>
#include <RAT/PosGen_Paint.hh>
#include <RAT/PosGen_Line.hh>
#include <RAT/PosGen_FillShell.hh>
#include <RAT/PosGen_Multipoint.hh>
#include <RAT/VertexGen_PhotonBomb.hh>
#include <RAT/VertexGen_LaserBall.hh>
#include <stdio.h>                   // for sprintf

#include <RAT/EventInfo.hh>

GLG4PrimaryGeneratorAction*
GLG4PrimaryGeneratorAction::fTheGLG4PrimaryGeneratorAction=0;

GLG4PrimaryGeneratorAction::
GLG4PrimaryGeneratorAction()
{
  if (fTheGLG4PrimaryGeneratorAction == 0) {
    fTheGLG4PrimaryGeneratorAction= this;
  }
  else {
    G4Exception(__FILE__, "Multiple Primary Generator Action", FatalException, "Error, more than one GLG4PrimaryGeneratorAction instantiated.\n"
    "Sorry, but this is a no-no because GLG4SteppingAction relies on\n"
    "GLG4PrimaryGeneratorAction::GetTheGLG4PrimaryGeneratorAction().\n"
    "This is yucky, I know -- please rewrite GLG4SteppingAction AND\n"
    "all main() programs so that constructor accepts a pointer to\n"
    "the GLG4PrimaryGeneratorAction you really want them to use.");
  }

  // initialize messenger and time fields
  fMessenger = new GLG4PrimaryGeneratorMessenger(this);
  fUniversalTime= 0.0;
  fUniversalTimeSincePriorEvent= 0.0;

  // initialize generator state
  GLG4param &db ( GLG4param::GetDB() );
  
  fEventWindow= db.GetWithDefault("gen.eventWindow",  0.*ns);

  // load up factories with known generators

  // Generic generators
  RAT::GlobalFactory<GLG4Gen>::Register("combo", 
					new RAT::Alloc<GLG4Gen, GLG4Gen_Combo>);

  // These generators are used by combo to make a "custom" generator
  // Rate generators
  RAT::GlobalFactory<GLG4TimeGen>::Register("uniform", 
					    new RAT::Alloc<GLG4TimeGen, 
					    GLG4TimeGen_Uniform>);
  RAT::GlobalFactory<GLG4TimeGen>::Register("poisson",
					    new RAT::Alloc<GLG4TimeGen, 
					    GLG4TimeGen_Poisson>);

  // Vertex generators
  RAT::GlobalFactory<GLG4VertexGen>::Register("gun", 
					      new RAT::Alloc<GLG4VertexGen, 
					      GLG4VertexGen_Gun>);
  RAT::GlobalFactory<GLG4VertexGen>::Register("gun2", 
					      new RAT::Alloc<GLG4VertexGen, 
					      GLG4VertexGen_Gun2>);
  RAT::GlobalFactory<GLG4VertexGen>::Register("gun3",
                                              new RAT::Alloc<GLG4VertexGen,
                                              GLG4VertexGen_Gun2>);
  RAT::GlobalFactory<GLG4VertexGen>::Register("pbomb", 
					      new RAT::Alloc<GLG4VertexGen, 
					      RAT::VertexGen_PhotonBomb>);
  RAT::GlobalFactory<GLG4VertexGen>::Register("lsbp", 
					      new RAT::Alloc<GLG4VertexGen, 
					      RAT::VertexGen_LaserBall>);

  // Position generators
  RAT::GlobalFactory<GLG4PosGen>::Register("point", new RAT::Alloc<GLG4PosGen, RAT::PosGen_Point>);
  RAT::GlobalFactory<GLG4PosGen>::Register("line",  new RAT::Alloc<GLG4PosGen, RAT::PosGen_Line>);
  RAT::GlobalFactory<GLG4PosGen>::Register("paint", new RAT::Alloc<GLG4PosGen, RAT::PosGen_Paint>);
  RAT::GlobalFactory<GLG4PosGen>::Register("fill",  new RAT::Alloc<GLG4PosGen, RAT::PosGen_Fill>);
  RAT::GlobalFactory<GLG4PosGen>::Register("fillshell", 
                                           new RAT::Alloc<GLG4PosGen, 
                                           RAT::PosGen_FillShell>);
  RAT::GlobalFactory<GLG4PosGen>::Register("multipoint",  
                                           new RAT::Alloc<GLG4PosGen,  
                                           RAT::PosGen_Multipoint>);

  fNeedReset = true;
}

GLG4PrimaryGeneratorAction::~GLG4PrimaryGeneratorAction()
{
}

void GLG4PrimaryGeneratorAction::SetEventWindow(double argEventWindow)
{
  fEventWindow= argEventWindow;
  
  GLG4param &db ( GLG4param::GetDB() );
  db["gen.eventWindow"]= argEventWindow;
}

void GLG4PrimaryGeneratorAction::AddGenerator(GLG4Gen *gen)
{
  fGenList.push(gen);
  fNeedReset = true;
}

// GeneratePrimaries (this is the interesting part!)
void GLG4PrimaryGeneratorAction::GeneratePrimaries(G4Event* argEvent)
{
  // Add RAT specific information here (FIXME in some alternate universe where
  // this program continues to be developed)
  argEvent->SetUserInformation(new RAT::EventInfo);

  
  if (fNeedReset) {
    // Need to reset the times of all the generators now that we've started
    // running.

    GLG4GenList temp;
    while (!fGenList.empty()) {
      GLG4Gen *gen = fGenList.top();
      fGenList.pop();
      gen->ResetTime();
      temp.push(gen);
    }
    
    fGenList = temp; // copy them all back
    fNeedReset = false;
  }

  // Find time increment until next event on queue
  // (priority queue will always have next event on top)
  if (fGenList.empty())
    G4Exception(__FILE__, "No Active Generator", FatalException, "GLG4PrimaryGeneratorAction: No generators selected!");
   
  double timeToNextEvent = fGenList.top()->NextTime();
  fUniversalTimeSincePriorEvent = timeToNextEvent;
  fUniversalTime += timeToNextEvent;

  // Offset time of all queued generators, so top generator's NextTime() == 0
  fGenList.SubtractTime(timeToNextEvent);

  // Add all events in event window, includes pileup and deferred tracks
  while (!fGenList.empty()
	 && fGenList.top()->NextTime() <= fEventWindow) {

    GLG4Gen *nextGen = fGenList.top();
    fGenList.pop();

    double vertexTime = nextGen->NextTime();
    nextGen->GenerateEvent(argEvent);

    if (nextGen->IsRepeatable()) {
      // Reset time relative to time to previous generator time
      // to get pilup correct
      nextGen->ResetTime(vertexTime);
      fGenList.push(nextGen);
    } else
      delete nextGen;
  }
}


void GLG4PrimaryGeneratorAction::
DeferTrackToLaterEvent(const G4Track * track)
{
  fGenList.push(new GLG4Gen_DeferTrack(track));
}
