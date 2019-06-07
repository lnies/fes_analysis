#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "aio_module_def.h"

 void timestring(char * string){
   time_t now;
   struct tm today;

   now = time(NULL);     // get date and time
   today = *localtime(&now);
   sprintf(string, "%02d.%02d.%d %02d:%02d:%02d", 
	  today.tm_mday, today.tm_mon +1, today.tm_year+1900,
          today.tm_hour, today.tm_min, today.tm_sec);
}

void binarystring(int value, char* string){
  string[0]='0';
  string[1]='b';
  for(int n=31;n>=0; n--){
      if(((value>>n)&0x1)==1) string[31-n+2]='1';
       else string[31-n+2]='0';
  }
  string[34]=0;
}


void AIO_what_to_name(int what, char *name){
  switch(what){
    case 11:     
      strcpy(name, "ntec     4ch BaF2");
    break;
    case 21:
      strcpy(name, "veto     8ch pl. Veto");
    break;
    case 30:
      strcpy(name, "v775N   16ch TDC");
    break;
    case 31:
      strcpy(name, "v775    32ch TDC");
    break;
    case 40:
      strcpy(name, "v785N   16ch PS ADC");
    break;
    case 41:
      strcpy(name, "v785    32ch PS ADC");
    break;
    case 50:
      strcpy(name, "v792N   16ch QDC");
    break;
    case 51:
      strcpy(name, "v792    32ch QDC");
    break;
    case 60:
      strcpy(name, "v965A    8ch DR QDC");
    break;
    case 61:
      strcpy(name, "v965    16ch DR QDC");
    break;
    case 89:
      strcpy(name, "IOL board, use this for LAM");
    break;
    case 90:
      strcpy(name, "sis3820 Sync module, use this for LAM");
    break;
    case 91:
      strcpy(name, "Bonn sync module, no LAM");
    break;
    case 110:
      strcpy(name, "v895    16ch LED or the CFD");
    break;
    case 160:
      strcpy(name, "sis3302  SIS3302 8ch 16bit FADC");
    break;
    case 170:
      strcpy(name, "Bonn Tracking Detectorsystem");
    break;
    case 180:
      strcpy(name, "AVM16 WIENER 16 ch ADC, 160MHz");
    break;
    default:
      strcpy(name, "unknown module");
    break;
  } 
}
