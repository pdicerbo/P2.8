#include "Vector.h"
#include "Random.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <mpi.h>

using namespace std;
using namespace PLMD;


class SimpleMD
{

int iv[32];
int iy;
int iset;
double gset;
bool write_positions_first;
bool write_statistics_first;
int write_statistics_last_time_reopened;
FILE* write_statistics_fp;

int modulo(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}


public:
SimpleMD(){
  for(int i=0;i<32;i++) iv[i]=0.0;
  iy=0;
  iset=0;
  gset=0.0;
  write_positions_first=true;
  write_statistics_first=true;
  write_statistics_last_time_reopened=0;
  write_statistics_fp=NULL;
}

private:

void 
read_input(FILE*   fp,
           double& temperature,
           double& tstep,
           double& friction,
           double& forcecutoff,
           double& listcutoff,
           int&    nstep,
           int&    nconfig,
           int&    nstat,
           bool&   wrapatoms,
           string& inputfile,
           string& outputfile,
           string& trajfile,
           string& statfile,
           int&    maxneighbours,
           int&    exchangestride,
           int&    idum)
{
  temperature=1.0;
  tstep=0.005;
  friction=0.0;
  forcecutoff=2.5;
  listcutoff=3.0;
  nstep=1;
  nconfig=10;
  nstat=1;
  maxneighbours=1000;
  idum=0;
  wrapatoms=false;
  statfile="";
  trajfile="";
  outputfile="";
  inputfile="";
  exchangestride=0;

  string line;

  line.resize(256);
  char buffer[256];
  char buffer1[256];

  while(fgets(buffer,256,fp)){
    line=buffer;
    for(int i=0;i<line.length();++i) if(line[i]=='#' || line[i]=='\n') line.erase(i);
    for(int i=line.length()-1;i>=0;--i){
      if(line[i]!=' ')break;
      line.erase(i);
    }
    buffer[0]=0;
    sscanf(line.c_str(),"%s",buffer);
    if(strlen(buffer)==0) continue;
    string keyword=buffer;
    if(keyword=="temperature")
      sscanf(line.c_str(),"%s %lf",buffer,&temperature);
    else if(keyword=="tstep")
      sscanf(line.c_str(),"%s %lf",buffer,&tstep);
    else if(keyword=="friction")
      sscanf(line.c_str(),"%s %lf",buffer,&friction);
    else if(keyword=="forcecutoff")
      sscanf(line.c_str(),"%s %lf",buffer,&forcecutoff);
    else if(keyword=="listcutoff")
      sscanf(line.c_str(),"%s %lf",buffer,&listcutoff);
    else if(keyword=="nstep")
      sscanf(line.c_str(),"%s %d",buffer,&nstep);
    else if(keyword=="nconfig")
    {
      sscanf(line.c_str(),"%s %d %s",buffer,&nconfig,buffer1);
      trajfile=buffer1;
    }
    else if(keyword=="nstat")
    {
      sscanf(line.c_str(),"%s %d %s",buffer,&nstat,buffer1);
      statfile=buffer1;
    }
    else if(keyword=="wrapatoms")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') wrapatoms=true;
    }
    else if(keyword=="maxneighbours")
      sscanf(line.c_str(),"%s %d",buffer,&maxneighbours);
    else if(keyword=="inputfile")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      inputfile=buffer1;
    }
    else if(keyword=="outputfile")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      outputfile=buffer1;
    }
    else if(keyword=="idum")
      sscanf(line.c_str(),"%s %d",buffer,&idum);
    else if(keyword=="exchangestride")
      sscanf(line.c_str(),"%s %d",buffer,&exchangestride);
    else{
      fprintf(stderr,"Unknown keywords :%s\n",keyword.c_str());
      exit(1);
    }
  }

  if(inputfile.length()==0){
      fprintf(stderr,"Specify input file\n");
      exit(1);
  }
  if(outputfile.length()==0){
      fprintf(stderr,"Specify output file\n");
      exit(1);
  }
  if(trajfile.length()==0){
      fprintf(stderr,"Specify traj file\n");
      exit(1);
  }
  if(statfile.length()==0){
      fprintf(stderr,"Specify stat file\n");
      exit(1);
  }
}

void read_natoms(const string & inputfile,int & natoms){
// read the number of atoms in file "input.xyz"
  FILE* fp=fopen(inputfile.c_str(),"r");
  fscanf(fp,"%d",&natoms);
  fclose(fp);
}

void read_positions(const string& inputfile,int natoms,vector<Vector>& positions,double cell[3]){
// read positions and cell from a file called inputfile
// natoms (input variable) and number of atoms in the file should be consistent
  FILE* fp=fopen(inputfile.c_str(),"r");
  char buffer[256];
  char atomname[256];
  fgets(buffer,256,fp);
  fscanf(fp,"%lf %lf %lf",&cell[0],&cell[1],&cell[2]);
  for(int i=0;i<natoms;i++){
    fscanf(fp,"%s %lf %lf %lf",atomname,&positions[i][0],&positions[i][1],&positions[i][2]);
// note: atomname is read but not used
  }
  fclose(fp);
}

void randomize_velocities(const int natoms,const double temperature,const vector<double>&masses,vector<Vector>& velocities,Random&random){
// randomize the velocities according to the temperature
  for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++)
      velocities[iatom][i]=sqrt(temperature/masses[iatom])*random.Gaussian();
}

void pbc(const double cell[3],const Vector & vin,Vector & vout){
// apply periodic boundary condition to a vector
  for(int i=0;i<3;i++){
    vout[i]=vin[i]-floor(vin[i]/cell[i]+0.5)*cell[i];
  }
}

void check_list(const int natoms,const vector<Vector>& positions,const vector<Vector>&positions0,const double listcutoff,
                const double forcecutoff,bool & recompute)
{
// check if the neighbour list have to be recomputed
  Vector displacement;  // displacement from positions0 to positions
  double delta2;        // square of the 'skin' thickness
  recompute=false;
  delta2=(0.5*(listcutoff-forcecutoff))*(0.5*(listcutoff-forcecutoff));
// if ANY atom moved more than half of the skin thickness, recompute is set to .true.
  for(int iatom=0;iatom<natoms;iatom++){
    for(int k=0;k<3;k++) displacement[k]=positions[iatom][k]-positions0[iatom][k];
    double s=0.0;
    for(int k=0;k<3;k++) s+=displacement[k]*displacement[k];
    if(s>delta2) recompute=true;
  }
}


void compute_list(const int natoms,const vector<Vector>& positions,const double cell[3],const double listcutoff,
                  vector<vector<int> >& list,MPI_Comm intracomm){
  int rank=0;
  int npe=1;
#ifdef SIMPLEMD_MPI
  MPI_Comm_rank(intracomm,&rank);
  MPI_Comm_size(intracomm,&npe);
#endif

  Vector distance;     // distance of the two atoms
  Vector distance_pbc; // minimum-image distance of the two atoms
  double listcutoff2;  // squared list cutoff
  listcutoff2=listcutoff*listcutoff;
  list.assign(natoms,vector<int>());

  Vector domain_side;
  int ndomains[3];
  for(int k=0;k<3;k++) ndomains[k]=int(cell[k]/listcutoff);
  std::vector<std::vector<int> > domains;
  domains.resize(ndomains[0]*ndomains[1]*ndomains[2]);
  for(int k=0;k<3;k++) domain_side[k]=cell[k]/ndomains[k];
  for(int iatom=0;iatom<natoms;iatom++){
    int idomain[3];
    for(int k=0;k<3;k++) idomain[k]=modulo(int((positions[iatom][k]/domain_side[k])),ndomains[k]);
    domains[(idomain[0]*ndomains[1]+idomain[1])*ndomains[2]+idomain[2]].push_back(iatom);
  }

  int count=0;
  for(int ix=0;ix<ndomains[0];ix++) for(int iy=0;iy<ndomains[1];iy++) for(int iz=0;iz<ndomains[2];iz++)
  if(count++%npe==rank )
  for(int jx=0;jx<ndomains[0];jx++) for(int jy=0;jy<ndomains[1];jy++) for(int jz=0;jz<ndomains[2];jz++)
  if( (ix-jx+ndomains[0])%ndomains[0]<=1 || (ix-jx+ndomains[0])%ndomains[0]>ndomains[0]-2)
  if( (iy-jy+ndomains[1])%ndomains[1]<=1 || (iy-jy+ndomains[1])%ndomains[1]>ndomains[1]-2)
  if( (iz-jz+ndomains[2])%ndomains[2]<=1 || (iz-jz+ndomains[2])%ndomains[2]>ndomains[2]-2)
  for(int iiatom=0;iiatom<domains[(ix*ndomains[1]+iy)*ndomains[2]+iz].size();iiatom++){
    int iatom=domains[(ix*ndomains[1]+iy)*ndomains[2]+iz][iiatom];
    for(int jjatom=0;jjatom<domains[(jx*ndomains[1]+jy)*ndomains[2]+jz].size();jjatom++){
      int jatom=domains[(jx*ndomains[1]+jy)*ndomains[2]+jz][jjatom];
      if(jatom<=iatom) continue;
      for(int k=0;k<3;k++) distance[k]=positions[iatom][k]-positions[jatom][k];
      pbc(cell,distance,distance_pbc);
// if the interparticle distance is larger than the cutoff, skip
      double d2=0; for(int k=0;k<3;k++) d2+=distance_pbc[k]*distance_pbc[k];
      if(d2>listcutoff2)continue;
      list[iatom].push_back(jatom);
    }
  }
}

void compute_forces(const int natoms,const vector<Vector>& positions,const double cell[3],
                    double forcecutoff,const vector<vector<int> >& list,vector<Vector>& forces,double & engconf,MPI_Comm intracomm)
{
  Vector distance;        // distance of the two atoms
  Vector distance_pbc;    // minimum-image distance of the two atoms
  double distance_pbc2;   // squared minimum-image distance
  double forcecutoff2;    // squared force cutoff
  Vector f;               // force
  double engcorrection;   // energy necessary shift the potential avoiding discontinuities

  forcecutoff2=forcecutoff*forcecutoff;
  engconf=0.0;
  for(int i=0;i<natoms;i++)for(int k=0;k<3;k++) forces[i][k]=0.0;
  engcorrection=4.0*(1.0/pow(forcecutoff2,6.0)-1.0/pow(forcecutoff2,3));
  for(int iatom=0;iatom<natoms-1;iatom++){
    for(int jlist=0;jlist<list[iatom].size();jlist++){
      int jatom=list[iatom][jlist];
      for(int k=0;k<3;k++) distance[k]=positions[iatom][k]-positions[jatom][k];
      pbc(cell,distance,distance_pbc);
      distance_pbc2=0.0; for(int k=0;k<3;k++) distance_pbc2+=distance_pbc[k]*distance_pbc[k];
// if the interparticle distance is larger than the cutoff, skip
      if(distance_pbc2>forcecutoff2) continue;
      double distance_pbc6=distance_pbc2*distance_pbc2*distance_pbc2;
      double distance_pbc8=distance_pbc6*distance_pbc2;
      double distance_pbc12=distance_pbc6*distance_pbc6;
      double distance_pbc14=distance_pbc12*distance_pbc2;
      engconf+=4.0*(1.0/distance_pbc12 - 1.0/distance_pbc6) - engcorrection;
      for(int k=0;k<3;k++) f[k]=2.0*distance_pbc[k]*4.0*(6.0/distance_pbc14-3.0/distance_pbc8);
// same force on the two atoms, with opposite sign:
      for(int k=0;k<3;k++) forces[iatom][k]+=f[k];
      for(int k=0;k<3;k++) forces[jatom][k]-=f[k];
    }
  }

#ifdef SIMPLEMD_MPI
  MPI_Allreduce(MPI_IN_PLACE,&forces[0][0],3*natoms,MPI_DOUBLE,MPI_SUM,intracomm);
  MPI_Allreduce(MPI_IN_PLACE,&engconf,1,MPI_DOUBLE,MPI_SUM,intracomm);
#endif

}

void compute_engkin(const int natoms,const vector<double>& masses,const vector<Vector>& velocities,double & engkin)
{
// calculate the kinetic energy from the velocities
  engkin=0.0;
  for(int iatom=0;iatom<natoms;iatom++)for(int k=0;k<3;k++){
    engkin+=0.5*masses[iatom]*velocities[iatom][k]*velocities[iatom][k];
  }
}


void thermostat(const int natoms,const vector<double>& masses,const double dt,const double friction,
                const double temperature,vector<Vector>& velocities,double & engint,Random & random){
// Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
  double c1,c2;
  c1=exp(-friction*dt);
  for(int iatom=0;iatom<natoms;iatom++){
    c2=sqrt((1.0-c1*c1)*temperature/masses[iatom]);
    for(int i=0;i<3;i++){
      engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
      engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
    }
  }
}

void write_positions(const string& trajfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
{
// write positions on file trajfile
// positions are appended at the end of the file
  Vector pos;
  FILE*fp;
  if(write_positions_first){
    fp=fopen(trajfile.c_str(),"w");
    write_positions_first=false;
  } else {
    fp=fopen(trajfile.c_str(),"a");
  }
  fprintf(fp,"%d\n",natoms);
  fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
  for(int iatom=0;iatom<natoms;iatom++){
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
    if(wrapatoms) pbc(cell,positions[iatom],pos);
    else for(int k=0;k<3;k++) pos[k]=positions[iatom][k];
    fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
  }
  fclose(fp);
}

void write_final_positions(const string& outputfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
{
// write positions on file outputfile
  Vector pos;
  FILE*fp;
  fp=fopen(outputfile.c_str(),"w");
  fprintf(fp,"%d\n",natoms);
  fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
  for(int iatom=0;iatom<natoms;iatom++){
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
    if(wrapatoms) pbc(cell,positions[iatom],pos);
    else for(int k=0;k<3;k++) pos[k]=positions[iatom][k];
    fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
  }
  fclose(fp);
}


void write_statistics(const string & statfile,const int istep,const double tstep,
                      const int natoms,const double engkin,const double engconf,const double engint){
// write statistics on file statfile
  if(write_statistics_first){
// first time this routine is called, open the file
    write_statistics_fp=fopen(statfile.c_str(),"w");
    write_statistics_first=false;
  }
  if(istep-write_statistics_last_time_reopened>100){
// every 100 steps, reopen the file to flush the buffer
    fclose(write_statistics_fp);
    write_statistics_fp=fopen(statfile.c_str(),"a");
    write_statistics_last_time_reopened=istep;
  }
  fprintf(write_statistics_fp,"%d %f %f %f %f %f\n",istep,istep*tstep,2.0*engkin/(3.0*natoms),engconf,engkin+engconf,engkin+engconf+engint);
}



public:
#ifdef SIMPLEMD_MPI
int main(FILE*in,FILE*out,MPI_Comm intercomm,MPI_Comm intracomm){
#else
int main(FILE*in,FILE*out){
#endif
  int            natoms;       // number of atoms
  vector<Vector> positions;    // atomic positions
  vector<Vector> velocities;   // velocities
  vector<double> masses;       // masses
  vector<Vector> forces;       // forces
  double         cell[3];      // cell size
  double         cell9[3][3];  // cell size

// neighbour list variables
// see Allen and Tildesey book for details
  vector< vector<int> >  list;         // neighbour list
  vector<Vector> positions0;   // reference atomic positions, i.e. positions when the neighbour list

// input parameters
// all of them have a reasonable default value, set in read_input()
  double      tstep;             // simulation timestep
  double      temperature;       // temperature
  double      friction;          // friction for Langevin dynamics (for NVE, use 0)
  double      listcutoff;        // cutoff for neighbour list
  double      forcecutoff;       // cutoff for forces
  int         nstep;             // number of steps
  int         nconfig;           // stride for output of configurations
  int         nstat;             // stride for output of statistics
  int         maxneighbour;      // maximum average number of neighbours per atom
  int         idum;              // seed
  bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
  string      inputfile;         // name of file with starting configuration (xyz)
  string      outputfile;        // name of file with final configuration (xyz)
  string      trajfile;          // name of the trajectory file (xyz)
  string      statfile;          // name of the file with statistics
  string      string;            // a string for parsing

  double engkin;                 // kinetic energy
  double engconf;                // configurational energy
  double engint;                 // integral for conserved energy in Langevin dynamics

  bool recompute_list;           // control if the neighbour list have to be recomputed

  int exchangestride;            // stride for replica exchange

  Random random;                 // random numbers stream

  read_input(in,temperature,tstep,friction,forcecutoff,
             listcutoff,nstep,nconfig,nstat,
             wrapatoms,inputfile,outputfile,trajfile,statfile,
             maxneighbour,exchangestride,idum);

#ifdef SIMPLEMD_MPI
  int rank;
  int npe;
  MPI_Comm_rank(intracomm,&rank);
  MPI_Comm_size(intracomm,&npe);
// avoid output in ranks>0
  if(rank>0){
    trajfile="/dev/null";
    statfile="/dev/null";
    outputfile="/dev/null";
  }
#endif



// number of atoms is read from file inputfile
  read_natoms(inputfile,natoms);

// write the parameters in output so they can be checked
  fprintf(stdout,"%s %s\n","Starting configuration           :",inputfile.c_str());
  fprintf(stdout,"%s %s\n","Final configuration              :",outputfile.c_str());
  fprintf(stdout,"%s %d\n","Number of atoms                  :",natoms);
  fprintf(stdout,"%s %f\n","Temperature                      :",temperature);
  fprintf(stdout,"%s %f\n","Time step                        :",tstep);
  fprintf(stdout,"%s %f\n","Friction                         :",friction);
  fprintf(stdout,"%s %f\n","Cutoff for forces                :",forcecutoff);
  fprintf(stdout,"%s %f\n","Cutoff for neighbour list        :",listcutoff);
  fprintf(stdout,"%s %d\n","Number of steps                  :",nstep);
  fprintf(stdout,"%s %d\n","Stride for trajectory            :",nconfig);
  fprintf(stdout,"%s %s\n","Trajectory file                  :",trajfile.c_str());
  fprintf(stdout,"%s %d\n","Stride for statistics            :",nstat);
  fprintf(stdout,"%s %s\n","Statistics file                  :",statfile.c_str());
  fprintf(stdout,"%s %d\n","Max average number of neighbours :",maxneighbour);
  fprintf(stdout,"%s %d\n","Seed                             :",idum);
  fprintf(stdout,"%s %d\n","Stride for replica exchange      :",exchangestride);
  fprintf(stdout,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));

// Setting the seed
  random.setSeed(idum);

// allocation of dynamical arrays
  positions.resize(natoms);
  positions0.resize(natoms);
  velocities.resize(natoms);
  forces.resize(natoms);
  masses.resize(natoms);
  list.resize(natoms);

// masses are hard-coded to 1
  for(unsigned i=0;i<natoms;++i) masses[i]=1.0;

// energy integral initialized to 0
  engint=0.0;

// positions are read from file inputfile
  read_positions(inputfile,natoms,positions,cell);

// velocities are randomized according to temperature
  randomize_velocities(natoms,temperature,masses,velocities,random);

// neighbour list are computed, and reference positions are saved
  compute_list(natoms,positions,cell,listcutoff,list,intracomm);

  int list_size=0;
  for(int i=0;i<list.size();i++) list_size+=list[i].size();
  fprintf(stdout,"List size: %d\n",list_size);
  for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];

// forces are computed before starting md
  compute_forces(natoms,positions,cell,forcecutoff,list,forces,engconf,intracomm);

// here is the main md loop
// Langevin thermostat is applied before and after a velocity-Verlet integrator
// the overall structure is:
//   thermostat
//   update velocities
//   update positions
//   (eventually recompute neighbour list)
//   compute forces
//   update velocities
//   thermostat
//   (eventually dump output informations)
  for(int istep=0;istep<nstep;istep++){
    thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

    for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
      velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];

    for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
      positions[iatom][k]+=velocities[iatom][k]*tstep;

// a check is performed to decide whether to recalculate the neighbour list
    check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list);
    if(recompute_list){
      compute_list(natoms,positions,cell,listcutoff,list,intracomm);
      for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];
      fprintf(stdout,"Neighbour list recomputed at step %d\n",istep);
      int list_size=0;
      for(int i=0;i<list.size();i++) list_size+=list[i].size();
      fprintf(stdout,"List size: %d\n",list_size);
    }

    compute_forces(natoms,positions,cell,forcecutoff,list,forces,engconf,intracomm);

    for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
      velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];

    thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

// kinetic energy is calculated
  compute_engkin(natoms,masses,velocities,engkin);

// eventually, write positions and statistics
    if((istep+1)%nconfig==0) write_positions(trajfile,natoms,positions,cell,wrapatoms);
    if((istep+1)%nstat==0)   write_statistics(statfile,istep+1,tstep,natoms,engkin,engconf,engint);

    if(exchangestride!=0 && (istep+1)%exchangestride==0){
      int nrep;
      int irep;
      MPI_Comm_size(intercomm,&nrep);
      MPI_Comm_rank(intercomm,&irep);
      double data[3];
      double partnerdata[3];
      data[0]=1.0/temperature;
      data[1]=engconf;
      data[2]=random.U01();
      int partner=irep+(((istep+1)/exchangestride+irep)%2)*2-1;
      if(partner<0) partner=0;
      if(partner>=nrep) partner=nrep-1;
      MPI_Status status;
      MPI_Sendrecv(data,3,MPI_DOUBLE,partner,123,partnerdata,3,MPI_DOUBLE,partner,123,intercomm,&status);

      double delta=(data[0]-partnerdata[0])*(data[1]-partnerdata[1]);
      double acceptance=exp(delta);
      if(acceptance>1.0) acceptance=1.0;
      fprintf(stderr,"### %d %d %d %lf %lf %lf %lf %lf\n", istep+1,irep,partner,data[0],data[1],partnerdata[0],partnerdata[1],acceptance);
// make sure both processes uses the same random number
      double r=0.0;
      if(partner<irep) r=data[2];
      else r=partnerdata[2];
      if(acceptance>r){
        vector<Vector> newpos(natoms);
        vector<Vector> newvel(natoms);
        MPI_Sendrecv(&positions[0],3*natoms,MPI_DOUBLE,partner,124,&newpos[0],3*natoms,MPI_DOUBLE,partner,124,intercomm,&status);
        MPI_Sendrecv(&velocities[0],3*natoms,MPI_DOUBLE,partner,125,&newvel[0],3*natoms,MPI_DOUBLE,partner,125,intercomm,&status);
        positions=newpos;
        velocities=newvel;
        double scale=sqrt(partnerdata[0]/data[0]);
        for(unsigned i=0;i<natoms;i++) for(unsigned k=0;k<3;k++) velocities[i][k]*=scale;
      }
      
    }

  }

  write_final_positions(outputfile,natoms,positions,cell,wrapatoms);

// close the statistic file if it was open:
  if(write_statistics_fp) fclose(write_statistics_fp);

  return 0;
}


};

int main(int argc,char*argv[]){
#ifdef SIMPLEMD_MPI
  MPI_Init(&argc,&argv);
  assert(argc>1);
  MPI_Comm intracomm;
  MPI_Comm intercomm;
  int nrep=argc-1;
  int npe,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&npe);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  assert(npe%nrep==0);
  MPI_Comm_split(MPI_COMM_WORLD,rank/(npe/nrep),rank,&intracomm);
  MPI_Comm_split(MPI_COMM_WORLD,rank%(npe/nrep),rank,&intercomm);
  int irep;
  MPI_Comm_rank(intercomm,&irep);
#endif
  SimpleMD smd;
  FILE* in=stdin;
#ifndef SIMPLEMD_MPI
  if(argc>1) in=fopen(argv[1],"r");
#endif
#ifdef SIMPLEMD_MPI
  in=fopen(argv[irep+1],"r");
  int r=smd.main(in,stdout,intercomm,intracomm);
#else
  int r=smd.main(in,stdout);
#endif
  if(argc>1) fclose(in);
#ifdef SIMPLEMD_MPI
  MPI_Finalize();
#endif
  return r;
 }


