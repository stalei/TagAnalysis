//  © Shahram Talei @ 2020 The University of Alabama
//you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation; either version 3 of the License, or
//(at your option) any later version.
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
// This code is distributed as is and there is no warranty or technical support
// © Shahram Talei @ 2019
// Reading stored data in a sage file
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <gsl/gsl_rng.h>
//#include <ctype.h>

////////////////////////////////////////////////////////////
//// Definitios and variables
char SageDir[500];
struct Path_Names
{
    char paths[100];
}*SageFilesPath;

int SageFilesCount;
int NumGalaxies;

// This structune holds all of the information from the Sage files
struct SageGalaxies
{
  int   Type;
  int   FileNr;
  long long   GalaxyIndex;
  int   HaloIndex;
  int   FOFHaloIndex;
  int   TreeIndex;

  int   SnapNum;
  int   CentralGal;
  float CentralMvir;

   //properties of subhalo at the last time this galaxy was a central galaaxy
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;

   //baryonic reservoirs
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;
  //metals
 float MetalsColdGas;
 float MetalsStellarMass;
 float MetalsBulgeMass;
 float MetalsHotGas;
 float MetalsEjectedMass;
 float MetalsICS;

  //misc
 float Sfr;
 float SfrBulge;
 float SfrICS;
 float DiskScaleRadius;
 float Cooling;
 float Heating;
 float LastMajorMerger;
 float OutflowRate;

 float infallMvir;  //infall properties
 float infallVvir;
 float infallVmax;
 float r_heat;

}*SageOutput;
//typedef struct SageGalaxies

////////////////////////////////////////////////////////////
//// Functions


int CountSageFiles(int snap)
{
         //int i=0;
         int n=0;
         char line[10000];
         char file1[10000];
         FILE *fd;
         sprintf(file1, "%s/sagefile_%03d.txt",SageDir, snap);
         fd = fopen(file1, "r");
         while(fgets(line, sizeof(line), fd)!=NULL)
         {
            n++;// 1;
         }
         fclose(fd);
         return n;
}
void ReadSageFNames(int snap,struct Path_Names *SageFile)
{
        //printf("Reading sage file-name(s) for snap:%d\n",snap);
    int i=0;
    char line[10000];
    char file1[10000];
    char ABSpath[10000];//2
    char model[500];
    FILE *fd;
    sprintf(file1, "%s/sagefile_%03d.txt", SageDir, snap);
    fd = fopen(file1, "r");
    while(EOF != fscanf(fd, "%s%*[^\n]", line))
    {
       strcpy(ABSpath,line);//2
       strcpy(model,strrchr(ABSpath,'/'));
       //strcpy(SageFilesPath[i].paths, line);//1
       if (model[0] == '/')
         memmove(model, model+1, strlen(model));
       sprintf(SageFile[i].paths,"%s/%s",SageDir,model);
       i++;
    }
   fclose(fd);
   return;
 }

 int ReadSageHeader(int FileCount,struct Path_Names *SageFile)
 {
      int Ntrees;
          int NtotGal;
          int totmal;
          int i;
          FILE *fd;
          char file1[1000];
          totmal = 0;
          for(i=0; i<FileCount; i++)
             {
                  sprintf(file1, "%s", SageFile[i].paths);
                  fd = fopen(file1, "rb");
                  if(NULL == fd)
                  {
                     printf("Cannot open sage file");
                     return(-1);
                  }
                  fread(&Ntrees, sizeof(int), 1, fd);
                  fread(&NtotGal, sizeof(int), 1, fd);
                  totmal = totmal + NtotGal;
                  fclose(fd);
             }
 return totmal;
 }

void LoadSageFiles(int snap)
{
        SageFilesCount = CountSageFiles(snap);
        SageFilesPath = (struct Path_Names*)malloc(SageFilesCount * sizeof(struct Path_Names));
        ReadSageFNames(snap,SageFilesPath);
        return;
    }

    void ReadSageModel(int FileCount,struct Path_Names *SageFile,struct SageGalaxies *Output)
    {
    int Ntrees;
               int NtotGal;
               int *galpertree;
               char file1[1000];
               int offset;
               int i;
               FILE *fd;
               offset = 0;
               for(i=0; i<FileCount; i++)
               {
                  sprintf(file1, "%s", SageFile[i].paths);
                  fd = fopen(file1, "rb");
                        if(NULL == fd)
                        {
                           printf("Cannot open sage for reading");
                           return;
                        }
                        fread(&Ntrees, sizeof(int), 1, fd);
                        fread(&NtotGal, sizeof(int), 1, fd);
                        galpertree = (int*)malloc(Ntrees*sizeof(int));
                        fread(&galpertree[0], sizeof(int), Ntrees, fd);
                        fread(&Output[offset], sizeof(struct SageGalaxies), NtotGal, fd);
                        fclose(fd);
                        offset = offset + NtotGal;
                }
                        free(galpertree);
                        return;
    }

void ReadSage(int snap)
{

printf("Extracting information from the sage file.\n");
LoadSageFiles(snap);
fflush(stdout);
NumGalaxies = ReadSageHeader(SageFilesCount,SageFilesPath);

if(( SageOutput = (struct SageGalaxies*)malloc(NumGalaxies * sizeof(struct SageGalaxies))) == NULL)
    {
        printf("Failed to allocate for target sage...");
        return;
    }
ReadSageModel(SageFilesCount, SageFilesPath, SageOutput);
return;
}


void PrintGalaxyInfo(struct SageGalaxies *Output,int index)
{
  int i;
  i=index;
printf("~~~~~~~~~~~~~~~~~~~~~~∮∮∮∮∮∮∮∮∮∮∮∮~~~~~~~~~~~~~~~~~~~~~~\n");
printf("Galaxy Information- galaxy:%d, Type:%d, FileNr:%d\n",i,Output[i].Type,Output[i].FileNr);
printf("GalIndex:%lld,HaloIndex:%d,TreeIndex:%d,CentralGal:%d,CentralMvir:%f\n",Output[i].GalaxyIndex,Output[i].HaloIndex,
Output[i].TreeIndex,Output[i].CentralGal,Output[i].CentralMvir);
printf("Pos(x,y,z):(%f,%f,%f)\nVel(vx,vy,vz):(%f,%f,%f)\nSpin:(%f,%f,%f),Len:%d\n",Output[i].Pos[0],Output[i].Pos[1],
Output[i].Pos[2],Output[i].Vel[0],Output[i].Vel[1],Output[i].Vel[2],Output[i].Spin[0],Output[i].Spin[1],Output[i].Spin[2],Output[i].Len);

printf("Mvir:%f,Rvir:%f,Vvir:%f,Vmax:%f,VelDisp:%f\n",Output[i].Mvir,Output[i].Rvir,Output[i].Vvir,Output[i].Vmax,Output[i].VelDisp);

printf("ColdGas:%f,StellarMass:%f,BulgeMass:%f,HotGas:%f\nEjectedMass:%f,BHMass:%f,ICS:%f\n",Output[i].ColdGas,Output[i].StellarMass,
Output[i].BulgeMass,Output[i].HotGas,Output[i].EjectedMass,Output[i].BlackHoleMass,Output[i].ICS);

printf("\nMetals:\nColdGas:%f,StellarMass:%f,BulgeMass:%f,HotGas:%f\nEjectedMass:%f,ICS:%f\n\n",Output[i].MetalsColdGas,
Output[i].MetalsStellarMass,Output[i].MetalsBulgeMass,Output[i].MetalsHotGas,Output[i].MetalsEjectedMass,Output[i].MetalsICS);
printf("Sfr:%f,SfrBulge:%f,SfrICS:%f\nRd:%f,Cooling:%f,Heating:%f\nLastMajorMerger:%f,OutflowRate:%f\n",Output[i].Sfr,
Output[i].SfrBulge,Output[i].SfrICS,Output[i].DiskScaleRadius,Output[i].Cooling,Output[i].Heating,Output[i].LastMajorMerger,Output[i].OutflowRate);

printf("infallMvir:%f,infallVvir:%f,infallVmax:%f,r_heat:%f\n",Output[i].infallMvir,Output[i].infallVvir,Output[i].infallVmax,Output[i].r_heat);
//printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
//printf("∮∮∮∮∮∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮∮∮∮∮∮\n");
printf("~~~~~~~~~~~~~~~~~~~~~~∮∮∮∮∮∮∮∮∮∮∮∮~~~~~~~~~~~~~~~~~~~~~~\n");

return;

}
void ExportGalaxy(struct SageGalaxies *Output, int id)
{
  printf("%g,%g,%g\n",Output[id].Pos[0], Output[id].Pos[1],Output[id].Pos[2]);
}

int main()
{
sprintf(SageDir,"/home/shahram/Desktop/Research/3_Tagging/TagAnalysis");
ReadSage(264);
int i;
for(i=0;i<NumGalaxies;i++)
  ExportGalaxy(SageOutput,i);

  return 0;
}
