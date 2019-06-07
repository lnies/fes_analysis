#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "tagger_lib.h"

int Taggerfile_class::initialized=0;

void Taggerfile_class::init(char * file){
  char version[16];

  in = fopen(file, "rb");	
  if(in==NULL){ 
    printf("%s doesn't exist!!!\n", file);
    exit(-2);
  }
  fread(version, 16, 1, in);
  if(strcmp(version,"Explora2PNDV0.1")!=0){
    printf("Fileversion incorrect. Should be \"Explora2PNDV0.1\", but is \"%s\"\n", version);
    exit(-2);
  }
  initialized=1;
}

int Taggerfile_class::read_one_event(unsigned int no_of_ev){
  if(initialized==1){
    for(int n=0; n<32; n++){
      beamphotons[n][0]=beamphotons[n][1]=cherenkovhits[n]=0;
    }
     if(0==fread((char *)&dataheader,sizeof(dataheader), 1, in)) return -13;
     if(!(dataheader.a == 0xa && dataheader.b == 0xb)){
       printf("In event %i: invalid dataheader; a: 0x%x, b: 0x%x\n", no_of_ev, dataheader.a, dataheader.b);
       return -1;
     }
     printf("Event: %i, CBEvent: %i, ProtoEvent: %i, #Beamphotons: %i, #Cherenkovhits: %i\n", 
           no_of_ev, dataheader.EventCB, dataheader.EventProto, dataheader.nTaggerhits, dataheader.nCherenkovhits);
     if(no_of_ev>dataheader.EventProto){
       printf("Number of event (%i) > NoE in taggerfile (%i), EXIT\n", no_of_ev, dataheader.EventProto);
       exit(-2);
     }
     if(no_of_ev<dataheader.EventProto){
       if(no_of_ev==0){
         printf("Number of event (%i) < NoE in taggerfile (%i). Synchronizing.\n", no_of_ev, dataheader.EventProto);
         return(-3);
       }
       else{
         printf("Number of event (%i) < NoE in taggerfile (%i). EXIT.\n", no_of_ev, dataheader.EventProto);
         exit(-3);
       }
     }
     if(dataheader.nTaggerhits>0){ 
       if(0==fread((char *)beamphotons,sizeof(double)*2*dataheader.nTaggerhits, 1, in)) return -13;
     }
     if(dataheader.nCherenkovhits>0){
       if(0==fread((char *)cherenkovhits,sizeof(double)*dataheader.nCherenkovhits, 1, in)) return -13;
     }
     return dataheader.nTaggerhits;
  }
  return -42;
}

double Taggerfile_class::get_taggertime(int number){
  if(initialized==0){
    printf("No Taggerfile given\n");
    return -1;
  }
  if(number>=0 && number<dataheader.nTaggerhits) return beamphotons[number][1];
   else{
       printf("taggertime number out of bounds! (number %i, range: %i-1)\n", number, dataheader.nTaggerhits);
    return -1;
  }
}

double Taggerfile_class::get_taggerenergy(int number){
  if(initialized==0){
    printf("No Taggerfile given\n");
    return -1;
  }
  if(number>=0 && number<dataheader.nTaggerhits) return beamphotons[number][0];
   else{
       printf("taggerenergy number out of bounds! (number %i, range: %i-1)\n", number, dataheader.nTaggerhits);
    return -1;
  }
}

double Taggerfile_class::get_cherenkovtime(int number){
  if(initialized==0){
    printf("No Taggerfile given\n");
    return -1;
  }
  if(number>=0 && number<dataheader.nCherenkovhits) return cherenkovhits[number];
   else{
       printf("Cherenkovtime number out of bounds! (number %i, range: %i-1)\n", number, dataheader.nCherenkovhits);
    return -1;
  }
}

int Taggerfile_class::get_no_taggerhits(void){
  if(initialized==0){
    return 0;
  }
  return dataheader.nTaggerhits;
}

int Taggerfile_class::get_no_cherenkovhits(void){
  if(initialized==0){
    return 0;
  }
  return dataheader.nCherenkovhits;
}
