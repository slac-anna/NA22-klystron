#include <iostream>
#include <fstream>
#include <string>

#include "generator.hh"

using namespace std;

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fParticleGun = new G4ParticleGun(1);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
   std::ifstream InFile;
   InFile.open("../NA22-3.txt", std::ifstream::in);

   int count = 0;
   
// data format (MKS): Z, R, Theta, bg_Z, bg_R, bg_Theta, Current, Time_of_Flight, Potential, Specie_Name, Particle_Number
   double Z, R, Theta, bg_Z, bg_R, bg_Theta, Current_Time, Potential, P_ID;
   std::string Specie;
   double X, Y, Z0;
   double mom_X, mom_Y, mom_Z, mom2;
   double mass = 9.1093837e-31;
   double q = 1.60e-19;
   double c = 299792458;
   double c2 = c * c ;
   double mc = mass * c;
   double mc2 = (mc * c)* (mc * c);
   double energy = 0;

//specie
   G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition *particle = particleTable->FindParticle("e-");

//read data
   while (!InFile.eof()) {

      InFile >> Z >> R >> Theta >> bg_Z >> bg_R >> bg_Theta >> Current_Time >> Potential >> Specie >>P_ID;

//position
      X = R * cos(Theta);
      Y = R * sin(Theta);
      Z0 = -(0.140335+0.0087)/2+(180.75-165.75)/1000;  //align with the structure 

      G4ThreeVector pos(X*m, Y*m, Z0*m);      

//direction
      mom_X = bg_R * cos(bg_Theta) * mc;
      mom_Y = bg_R * sin(bg_Theta) * mc;
      mom_Z = bg_Z * mc;
 
      G4ThreeVector mom(mom_X, mom_Y, mom_Z);

//kinetic energy
      mom2 = mom_X * mom_X + mom_Y * mom_Y + mom_Z * mom_Z ;
      energy = ( c * sqrt(mass * mass * c2 + mom2) - mass * c2 ) / q;
      
      fParticleGun->SetParticlePosition(pos);
      fParticleGun->SetParticleMomentumDirection(mom);
      fParticleGun->SetParticleEnergy(energy*eV);    
      fParticleGun->SetParticleDefinition(particle);

    count++;
    
    std::cout<< " index " << P_ID << " pos " << pos << " mom " << mom  << " energy " << energy << std::endl;

    fParticleGun->GeneratePrimaryVertex(anEvent);  //an event having 720 particles?
   }

   std::cout<< " Total " << count << " particles generated " <<std::endl;

}





