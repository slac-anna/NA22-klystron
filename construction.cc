#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{

/////////////////////////////materials
    G4NistManager *nist = G4NistManager::Instance();

    G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material *CuMat   = nist->FindOrBuildMaterial("G4_Cu");
    G4Material *SSMat   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material *PMMat   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material *VaMat   = nist->FindOrBuildMaterial("G4_Galactic");  

    //define vacuum using a gas at low density or predefined Galactic 
    G4double atomicNumber = 1.; 
    G4double massOfMole = 1.008*g/mole; 
    G4double density = 1.e-25*g/cm3; 
    G4double temperature = 2.73*kelvin; 
    G4double pressure = 3.e-18*pascal; 
    //G4Material *Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure); 
    

/////////////////////////////dimensions
    G4double cr_l         = 2*3.1015*mm;   //cavity regular cell length
    G4double cr_r         = 6.916420*mm;   //cavity regular cell radius
    G4double cr_g         = 2*1.05*mm;     //cavity gap from the disk tip to disk tip
    G4double rounding1    = 0.127*mm;
    G4double disk_r       = 0.70*mm;       //nose radius
    G4double out_body_r   = 7.9375*mm;     //cavity body outer radius

    G4double bp_r         = 2.1*mm;        //beampipe inner radius
    G4double bp_up_z      = 8.7*mm;        //cavity left side beampipe z upto the center of nose4-5
    G4double bp_down_z    = 13.75*mm;      //cavity right side beampipe z at the nose5 tailpipe straight end

    G4double pipe_taper_mid_z  = 24.765*mm;   //tailpipe end z not including rounding tip
    G4double pipe_taper_mid_r  = 3.922146*mm; //nose5 tailpipe end flat surface outer radius

    G4double gap_z1 = 29.419909*mm;            
    G4double gap_z2 = 29.845*mm;            //gap    
    G4double gap_z3 = 29.337*mm;            //gap
    G4double gap_z4 = 28.575*mm;            //gap   
    G4double gap_z5 = 36.195*mm;            //gap  
    G4double gap_z6 = 31.877*mm;            //gap 
    G4double gap_z7 = 31.115*mm;            //gap                   
    G4double gap_z8 = 42.545*mm;            //gap  

    G4double gap_r1 = 4.69218*mm;                                   
    G4double gap_r2 = 5.193368*mm;          //gap                        
    G4double gap_r3 = 5.842*mm;             //gap                    
    G4double gap_r4 = 6.35*mm;              //gap            
    G4double gap_r5 = 12.7*mm;              //gap
    G4double gap_r6 = 10.668*mm;            //gap            
    G4double gap_r7 = 9.906*mm;   
    G4double gap_r8 = 8.255*mm;    
    G4double gap_r9 = 7.493*mm;   
    G4double gap_out_r=15.875*mm;

    G4double collector_end_z=120.015*mm;
    G4double collector_r = gap_r8;
    G4double collector_out_r = 12.7*mm;

    //simulation domain:
    G4double xWorld = 20.*mm;
    G4double yWorld = 20.*mm;
    G4double zWorld = 80.*mm;


///////////////////////////////////structure starting point, change particle data z value
//structure starting at z=-(140.335+8.7)/2*mm in Geant4 model correcponding to the position z= 165.75mm in drawing based on the z=0 at the 1st cavity center
    G4double offset1 = 8.7*mm; 
    G4double offset = offset1 - (140.335+8.7)*mm/2;  

    G4Box *solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0,G4ThreeVector(0., 0., 0), logicWorld, "physWorld", 0, false, 0, true);


////////////////////////////solid body****************************************************************************** 
    //part 1: nose 4-5, cavity cell, nose5 tailpile 
    //cavity cell block
    G4int numRZ1 = 9;
    G4double z1 = -bp_up_z, r1 = bp_r;
    G4double z2 = -cr_l/2, r2 = r1;
    G4double z3 = z2, r3 = cr_r;
    G4double z4 = -z3, r4 = r3; 
    G4double z5 = -z2, r5 = r2;
    G4double z6 = bp_down_z, r6 = r5;
    G4double z7 = pipe_taper_mid_z, r7 = pipe_taper_mid_r;
    G4double z8 = z7, r8 = out_body_r;
    G4double z9 = z1, r9 = r8;

    G4double rc[] = {r1, r2, r3, r4, r5, r6, r7, r8, r9};
    G4double zc[] = {z1, z2, z3, z4, z5, z6, z7, z8, z9};
    
    G4VSolid* solidCavity = new G4GenericPolycone("aPloycone1", 0.*deg, 360.*deg, numRZ1, rc, zc);
    G4LogicalVolume *logicCavity = new G4LogicalVolume(solidCavity, CuMat, "logicalCavity");
    G4VPhysicalVolume *phyCavity = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicCavity, "physCavity", logicWorld, false, 0, true);

    //cavity left nose cone
    G4double disk_r0    = bp_r + disk_r;         //origin r
    G4double disk_z0    = -(cr_g/2 + disk_r);    //origin z

    G4int numd = 9;
    G4int numnc = numd + 2;
    G4double rnc[numnc];
    G4double znc[numnc];

    G4double pi = atan(1.)*4.;
    G4double dth = pi/(numd - 1);
    for (G4int i = 0; i < numd; i++)
      {
	G4double angle = i * dth;
	rnc[i] = disk_r0 - disk_r*cos(angle);
	znc[i] = disk_z0 + disk_r*sin(angle);
      }

    rnc[numd]   = rnc[numd-1]; znc[numd]   = z2;
    rnc[numd+1] = bp_r;        znc[numd+1] = z2;

    G4VSolid* solidLeftNC = new G4GenericPolycone("aPloycone2", 0.*deg, 360.*deg, numnc, rnc, znc);
    G4LogicalVolume *logicLeftNC = new G4LogicalVolume(solidLeftNC, CuMat, "logicalLeftNC");
    G4VPhysicalVolume *phyLeftNC = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicLeftNC, "physLeftNC", logicWorld, false, 0, true);

    //cavity right nose cone
    disk_z0    = cr_g/2+disk_r;      //origin z

    G4double rncr[numnc];
    G4double zncr[numnc];

    for (G4int i = 0; i < numd; i++)
      {
	G4double angle = i * dth;
	rncr[i] = disk_r0 - disk_r*cos(angle);
	zncr[i] = disk_z0 - disk_r*sin(angle);
      }

    rncr[numd]   = rncr[numd-1]; zncr[numd]   = z5;
    rncr[numd+1] = bp_r;        zncr[numd+1] = z5;

    G4VSolid* solidRightNC = new G4GenericPolycone("aPloycone3", 0.*deg, 360.*deg, numnc, rncr, zncr);
    G4LogicalVolume *logicRightNC = new G4LogicalVolume(solidRightNC, CuMat, "logicalRightNC");
    G4VPhysicalVolume *phyRightNC = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicRightNC, "physRightNC", logicWorld, false, 0, true);


//******************************************************************************************************
   //Part 2: gap between nose5 tailpipe and collector, collector
   //gap SS cavity stack
    G4int numRZ2 = 16;
    G4double gz0 = z7, gr0 = r7;
    G4double gz1 = gap_z1, gr1 = gap_r1;
    G4double gz2 = gap_z2, gr2 = gap_r2;
    G4double gz3 = gz2, gr3 = gap_r3; 
    G4double gz4 = gap_z3, gr4 = gap_r4;
    G4double gz5 = gap_z4, gr5 = gr4;
    G4double gz6 = gz5, gr6 = gap_r5;
    G4double gz7 = gap_z5, gr7 = gr6;
    G4double gz8 = gz7, gr8 = gap_r6;
    G4double gz9 = gap_z6, gr9 = gr8;
    G4double gz10 = gap_z7, gr10 = gap_r7;
    G4double gz11 = gz10, gr11 = gap_r8;
    G4double gz12 = gz9, gr12 = gap_r9;
    G4double gz13 = gap_z8, gr13 = gr12;
    G4double gz14 = gz13, gr14 = gap_out_r;
    G4double gz15 = gz0, gr15 = gr14;

    G4double grc[] = {gr0, gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9, gr10, gr11, gr12, gr13, gr14, gr15};
    G4double gzc[] = {gz0, gz1, gz2, gz3, gz4, gz5, gz6, gz7, gz8, gz9, gz10, gz11, gz12, gz13, gz14, gz15};

    G4VSolid* solidGap = new G4GenericPolycone("aPloycone4", 0.*deg, 360.*deg, numRZ2,grc, gzc);
    G4LogicalVolume *logicGap = new G4LogicalVolume(solidGap, SSMat, "logicalGap");
    G4VPhysicalVolume *phyGap = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicGap, "physGap", logicWorld, false, 0, true);

    //collector
    G4int numRZ3 = 7;
    G4double crc[] = {gap_r9, gap_r9, 0, 0, 6.35*mm, collector_out_r, collector_out_r};
    G4double czc[] = {gap_z8, 120.015*mm, 131.723257*mm, 140.335*mm, 140.335*mm, 133.985*mm, gap_z8};

    G4VSolid* solidCollector = new G4GenericPolycone("aPloycone5", 0.*deg, 360.*deg, numRZ3, crc, czc);
    G4LogicalVolume *logicCollector = new G4LogicalVolume(solidCollector, CuMat, "logicalCollector");
    G4VPhysicalVolume *phyCollector = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicCollector, "physCollector", logicWorld, false, 0, true);


//*********************************************************************************************************
    //part3: vacuum
    G4int numRZt = 42; 

    G4double trc[] = {0, r1, rnc[10], rnc[0], rnc[1], rnc[2], rnc[3], rnc[4], rnc[5], rnc[6], rnc[7], rnc[8], rnc[9], r3, r4, rncr[9], rncr[8], rncr[7], rncr[6], rncr[5], rncr[4], rncr[3], rncr[2], rncr[1], rncr[0], rncr[10], r6, r7, gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9, gr10, gr11, gr12, gr12, 0};
    G4double tzc[] = {z1, z1, znc[10], znc[0], znc[1], znc[2], znc[3], znc[4], znc[5], znc[6], znc[7], znc[8], znc[9], z3, z4, zncr[9], zncr[8], zncr[7], zncr[6], zncr[5], zncr[4], zncr[3], zncr[2], zncr[1], zncr[0], zncr[10], z6, z7, gz1, gz2, gz3, gz4, gz5, gz6, gz7, gz8, gz9, gz10, gz11, gz12, 120.015*mm, 131.723257*mm};

    G4VSolid* solidVacuum = new G4GenericPolycone("aPloycone6", 0.*deg, 360.*deg, numRZt, trc, tzc);
    G4LogicalVolume *logicVacuum = new G4LogicalVolume(solidVacuum, VaMat, "logicalVacuum");
    G4VPhysicalVolume *phyVacuum = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicVacuum, "physVacuum", logicWorld, false, 0, true);

//******************************************************************************************************
   //part4: poles
    G4double porc[] = {7.9375*mm, 7.9375*mm, 16.76*mm};
    G4double poz1[] = {-5.753*mm, -3.753*mm, -3.753*mm};

    G4VSolid* solidPoleLeft = new G4GenericPolycone("aPloycone7", 0.*deg, 360.*deg, 3, porc, poz1);
    G4LogicalVolume *logicPoleLeft = new G4LogicalVolume(solidPoleLeft, SSMat, "logicPoleLeft");
    for(G4int j=0; j<2; j++)
     {
       G4VPhysicalVolume *phyPoleLeft = new G4PVPlacement(0, G4ThreeVector(0., 0., (offset+j*17.5*mm)), logicPoleLeft, "physPoleLeft", logicWorld, false, j, true);
     }

    G4double poz2[]= {5.753*mm, 3.753*mm, 3.753*mm};
    G4VSolid* solidPoleRight = new G4GenericPolycone("aPloycone8", 0.*deg, 360.*deg, 3, porc, poz2);
    G4LogicalVolume *logicPoleRight = new G4LogicalVolume(solidPoleRight, SSMat, "logicPoleRight");
    for(G4int j=0; j<2; j++)
     {
       G4VPhysicalVolume *phyPoleRight = new G4PVPlacement(0, G4ThreeVector(0., 0., (offset+j*10.*mm)), logicPoleRight, "physPoleRight", logicWorld, false, j, true);
     }

    G4double poz3[] = {21.753*mm, 23.753*mm, 23.753*mm};
    G4VSolid* solidPoleLeft3 = new G4GenericPolycone("aPloycone9", 0.*deg, 360.*deg, 3, porc, poz3);
    G4LogicalVolume *logicPoleLeft3 = new G4LogicalVolume(solidPoleLeft3, SSMat, "logicPoleLeft3");
    G4VPhysicalVolume *phyPoleLeft3 = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicPoleLeft3, "physPoleLeft3", logicWorld, false, 0, true);

//******************************************************************************************************

//part 5: PM magnets
    G4double magr1[] = {7.9375*mm, 7.9375*mm, 18*mm, 18*mm};
    G4double magz1[] = {-8.7*mm, -5.753*mm, -3.753*mm, -8.7*mm};
    G4VSolid* solidMag1Left = new G4GenericPolycone("aPloycone9", 0.*deg, 360.*deg, 4, magr1, magz1);
    G4LogicalVolume *logicMag1Left = new G4LogicalVolume(solidMag1Left, PMMat, "logicMag1Left");
    G4VPhysicalVolume *phyMag1Left = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicMag1Left, "physMag1Left", logicWorld, false, 0, true);

    G4double magr2[] = {16.76*mm, 7.9375*mm, 7.9375*mm, 16.76*mm};
    G4double magz2[]= {3.753*mm, 5.753*mm, 11.753*mm, 13.753*mm};
    G4VSolid* solidMag2Right = new G4GenericPolycone("aPloycone10", 0.*deg, 360.*deg, 4, magr2, magz2);
    G4LogicalVolume *logicMag2Right = new G4LogicalVolume(solidMag2Right, PMMat, "logicMag2Right");
    for(G4int j=0; j<2; j++)
     {
       G4VPhysicalVolume *phyMag2Right = new G4PVPlacement(0, G4ThreeVector(0., 0., (offset+j*10.*mm)), logicMag2Right, "physMag2Right", logicWorld, false, j, true);
     }

//***********************************************************************************************************
    //part 6: detectors
    //detector above the cavity cell body for radiation monitor
    G4Tubs *solidDetector = new G4Tubs("solidDetector", out_body_r, (out_body_r+2), 3.753, 0., 360.*deg);
    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicDetector");
    G4VPhysicalVolume *physDetector = new G4PVPlacement(0, G4ThreeVector(0., 0., offset), logicDetector, "physDetector", logicWorld, false, 0, true);

    //fScoringVolume = logicCavity;

    fScoringVolume = logicCollector;
    
    return physWorld;

}

void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new
    MySensitiveDetector("SensitiveDetector");

//add B-field
/*
    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0., 0., 1.*tesla));

    G4FieldManager* fieldManager = new G4FieldManager();

    fieldManager->SetDetectorField(magField);

    fieldManager->CreateChordFinder(magField);  //use default method

    G4bool pushToContained = true;

    logicDetector -> SetFieldManager(fieldManager, pushToContained);

    //G4AutoDelete::Register(magField);

    //G4AutoDelete::Register(fieldManager);
*/
    logicDetector->SetSensitiveDetector(sensDet);
}
