#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "aw_lib.h"
#include "aio_module_def.h"


int ReadSystem_class::init_readout(char * file, int ext_realign_first_noe){
  char *pointer;

  realign_first_noe=ext_realign_first_noe;

  if(strstr(file,".dat")!=NULL){
    strncpy(filename, file, 99);
    single_file=1;
    if(realign_first_noe==42){ // ok, multiple files starting with that file
      realign_first_noe=1; // it must be realigned
      single_file=0;       // multiple files
      pointer=strstr(file,".dat");  // find  .dat
      *pointer=0;                   // cut off .dat
      pointer= strrchr(file,'.');  // find last .
      file_number=atoi(pointer+1)-1;  // digits after . are our file_number
      *pointer=0;                   // new filename without file_number
      strncpy(filename, file, 99);
    }
  }
  else{  // 
    strncpy(filename, file, 99);
    single_file=0;
    file_number=-1;
  }
  old_no_of_boards=-1;
  if(open_file()!=1){
    printf("Error with opening file!\n");
    exit(0);
  }
  for(int n=0; n<no_of_boards; n++){
    printf("Initializing module %i\n", n);
    switch(what[n]){
      case NTEC_MULTI:
      case VETO_MULTI:
      case V775N_TDC:
      case V775_TDC:
      case V785N_PSADC:
      case V785_PSADC:
      case V792N_QDC:
      case V792_QDC:
      case V965N_DRQDC:
      case V965_DRQDC:
      case IOL:
        read_modules[n] = new ReadRaw_class();
        break;
      case V895_LED:
        read_modules[n] = new ReadRaw_class_BLANK();
        break;
      case SIS3820_SYNC:
        read_modules[n] = new ReadRaw_class_SIS3820();
        break;
      case BONN_SYNC:
        read_modules[n] = new ReadRaw_class_BONNSYNC();
        break;
      case BONN_TRACK:
        read_modules[n] = new ReadRaw_class_BONNTRACK();
        break;
      case SIS3302_FADC:
        read_modules[n] = new ReadRaw_class_SIS3302();
        break;
      case AVM16_FADC:
        read_modules[n] = new ReadRaw_class_AVM16();
        break;
      default:
	printf("Unknown module: %i!\n", what[n]);
	exit (0);        
    }
    
    read_modules[n]->init(n, in, what[n], realign_first_noe);
  }   
} // ReadSystem_class::init_readout

int ReadSystem_class::open_file(void){
  char ververgl[200];
  int n, OK;
  
  if(single_file==0){ // different files
    file_number++;
    sprintf(file, "%s.%03i.dat", filename, file_number);
  }
  else{
    file_number=1;
    sprintf(file, "%s", filename);    
  }
  printf("\nReading in file %s\n", file);
  in = fopen(file, "rb");	
  if(in==NULL){ 
    printf("%s doesn't exist!!!\n", file);
    return(-2);
  }
  fseek(in, 0, 2); //setze Zeiger auf Ende der Datei
  no_of_int=ftell(in)/4; // Dateilaenge/4 pro header, data word, trail

  fseek(in, 0, 0); //setze Zeiger auf Anfang der Datei

  n=0;
  do{
    fread(&ververgl[n], 1, 1, in);
    n++;
  }while(ververgl[n-1]!=0);
  printf("The version of the readout was: %s (expected: %s)\n\n", 
            ververgl, AIO_versionsstr);
  n=0;
  do{
    fread(&gen_starttime[n], 1, 1, in);
    n++;
  }while(gen_starttime[n-1]!=0);
  printf("The readout was started: %s\n", 
         gen_starttime);
  n=0;
  do{
    fread(&starttime[n], 1, 1, in);
    n++;
  }while(starttime[n-1]!=0);
  printf("The file was started: %s\n", 
         starttime);

  n=0;
  do{
    fread(&stoptime[n], 1, 1, in);
    n++;
  }while(stoptime[n-1]!=0);
  printf("The readout was stopped: %s\n\n", 
          stoptime);

  fread(&no_of_boards, sizeof(no_of_boards), 1, in);
  printf("%i boards were read out.\n", no_of_boards);
  while(no_of_boards==0);
  fread(what, sizeof(what), 1, in); 
  fread(geos, sizeof(geos), 1, in); 
  for(int n=0; n<no_of_boards; n++){
    AIO_what_to_name(what[n], ververgl);
    printf("Board no %i is a %s (%i, geo %i)\n", 
	   n, ververgl, what[n], geos[n]);
    if(geos[n]>43) geos[n]=43;
    geo_to_bnr[ geos[n] ]=n;
  }
  OK=1;
  if(old_no_of_boards==-1){  // first time
    old_no_of_boards=no_of_boards;
    for(n=0; n<16; n++){
      old_what[n]=what[n];
      old_geos[n]=geos[n];
    }
  }
  else{  // same values as last time?
    if(old_no_of_boards!=no_of_boards){
      OK=0;
      printf("Error! Number of boards has changed between files!\n");
    }
    for(n=0; n<16; n++){
      if(old_what[n]!=what[n]){
        OK=0;
        printf("Error! Kind of board has changed between files!\n");
      }
      if(old_geos[n]!=geos[n]){
        OK=0;
        printf("Error! Geo address has changed between files!\n");
      }
    }
  }
  return OK;
}  // ReadSystem_class::open_file

unsigned int ReadSystem_class::get_value(
                                 int board, int channel, int event){
  //  printf("int: %i\n", channel);
  return read_modules[board]->getval(channel, event);
} // ReadSystem_class::get_value

unsigned int ReadSystem_class::get_trace(
                                 int board, int channel, unsigned int *array){
  return read_modules[board]->get_trace(channel, array);
} // ReadSystem_class::get_trace

int ReadSystem_class::get_multihitno(
                                 int module, int channel){
  return read_modules[module]->get_multihitno(channel);
} // ReadSystem_class::get_multihitno

int ReadSystem_class::get_further_parameter(
                                 int module, int number){
  return read_modules[module]->get_further_parameter(number);
} // ReadSystem_class::get_further_parameter
int ReadSystem_class::get_raw_ch_enabled(
                                 int module, int channel){
  return read_modules[module]->get_raw_ch_enabled(channel);
} // ReadSystem_class::get_raw_ch_enabled


int ReadSystem_class::read_one_event(unsigned int & no_of_ev){
    int m;
  for(int n=0; n<no_of_boards; n++){
//  printf("System reading out module %i\n", n);
    m= read_modules[n]->read_one_event(no_of_int, file_number, no_of_ev);
    if(m<0){ // realign noe
      if(realign_first_noe==1){
        printf("EXT: Realigning NOE to %i\n", -m);
        if(file_number!=0) m*=-1;
        printf("%i is realigned with offset %i\n", no_of_ev, m);
        no_of_ev+=(m*-1);
        printf("New no_of_ev: %i\n", no_of_ev, m);
        m=1;
	realign_first_noe=0;
      }
      else{ // do relignment only once
	printf("Realignment was done before, now it is an error!\n");
        m=2;
      }
    }     
    if(m==0 && single_file==0){ // End of file & multiple files
      printf("End of file %s\n", file);
      fclose(in);  // close old file
      if(open_file()!=1){ // open next file
	return m;         // if not ok
      }
       else{               // OK
        for(int x=0; x<no_of_boards; x++){
          read_modules[x]->set_filepointer(in);
          read_modules[x]->read_further_parameters();
	}
        n=0;
//        no_of_int=0;
        m= read_modules[n]->read_one_event(no_of_int, file_number, no_of_ev);
        if(m<0){ // realign noe
          if(realign_first_noe==1){
            printf("EXT: Realigning NOE to %i\n", -m);
            if(file_number!=0) m*=-1;
            printf("%i is realigned with offset %i\n", no_of_ev, m);
            no_of_ev+=(m*-1);
            printf("New no_of_ev: %i\n", no_of_ev, m);
            m=1;
	    realign_first_noe=0;
          }
           else{ // do relignment only once
	    printf("Realignment was done before, now it is an error!\n");
            m=2;
          }
        }     
      }
    }
//    printf("System read out module %i, return_value= %i\n", n, m);
    if(m==0) break;
    if(m==2) break;
  }
  return m;
} // ReadSystem_class::read_one_event

unsigned int ReadSystem_class::get_noe(int module){
  return  read_modules[module]->get_noe();
}

int ReadSystem_class::get_maxchannels(int module){
  return  read_modules[module]->get_maxchannels();
}

int ReadSystem_class::get_maxevents(int module){
  return  read_modules[module]->get_maxevents();
}


unsigned int ReadRaw_class::wie_oft=1000000;
unsigned int ReadRaw_class::events=0;


unsigned int ReadRaw_class::getval(int channel, int event) const{
  if(channel>=0 && channel<maxchannels){
    //    printf("TT %i: %i\n", channel ,data[channel] );
    return(data[channel]);
  }
   else{
    printf("Channels out of range!\n");
    return 0;
  }
} //  ReadRaw::getval

int ReadRaw_class::get_multihitno(int channel){
   return 1;
} //  ReadRaw::get_multihit_no

int ReadRaw_class::get_raw_ch_enabled(int channel){
   return 0;
} //  ReadRaw::get_raw_ch_enabled

int ReadRaw_class::get_further_parameter(int number){
  if(number<=(further_parameters[0]&0xff))
    return further_parameters[number];
  else{
      printf("Module %i hasn't got a parameter no %i\n", 
        module_number, number);
      exit(0);
  }
} //  ReadRaw::get_further_parameter

unsigned int ReadRaw_class::get_trace(int channel, unsigned int *array) const{
  if(channel>=0 && channel<maxchannels){
    array[0]=data[channel];
    return 0;
  }
   else{
    printf("Channels out of range!\n");
    return(1);
  }  
}

void ReadRaw_class::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=1;
  if(what%2==1){ // is equal to a large module = 32 channel
    maxchannels=32;
    data = new unsigned int[32];
  }
  else{ // 16channel
    maxchannels=16;
    data = new unsigned int[16];
  }
   read_further_parameters();

}  // ReadRaw_class::init

void ReadRaw_class::set_filepointer(FILE * ext_in){
  in=ext_in;
}

void ReadRaw_class::read_further_parameters(void){
  int no_of_par=0;
  int ext_module_number=0;

  fread(&further_parameters[0], sizeof(further_parameters[0]), 1, in); 

  printf("Parameter 0: 0x%x\n", further_parameters[0]);
  no_of_par=further_parameters[0]&0xff;
  ext_module_number=((further_parameters[0]>>8)&0xff);
  if(ext_module_number!=module_number){
      printf("Wrong module number %i (expected %i, raw: 0x%x)\n", 
              ext_module_number, module_number, further_parameters[0]);   
    exit(0);
  }
  if(no_of_par>19){
    printf("too many extended parameters (0x%x) in board %i!\n", 
    further_parameters[0], module_number);
    exit(0);
  }
  for(int n=1; n<=no_of_par; n++){
    fread(&further_parameters[n], sizeof(further_parameters[0]), 1, in); 
    printf("Parameter %i: %i (0x%x)\n", n, 
              further_parameters[n], further_parameters[n]);
  }
  events=0;  // new file, events per file have to be reset;
} // ReadRaw_class::read_further_parameters


int ReadRaw_class::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word
  unsigned int dec_geo, dec_tow, dec_noch, dec_data, dec_ch; 
                            // decoded event informations
// printf("Reading out normal module %i (maxch %i)\n", module_number, maxchannels);
  memset(data, 0, sizeof(data[0])*maxchannels);    // set all to 0 
  // Header  **********************
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  dec_tow=((dataword>>24) & 0x7);  // decode type of dataword
  if(dec_tow!=2){ // not a header
    printf("H: %i (0x%x)\n", dec_tow, dataword);
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0);
    dec_tow=((dataword>>24) & 0x7);  // decode type of dataword
    printf("After H: %i (0x%x)\n", dec_tow, dataword);
    exit(0);
  }
  dec_geo = (dataword>>27); // decode board number
  dec_noch = ((dataword>>8) & 0x3f);     // decode no of events in board
  if(dec_noch>maxchannels){
    printf("ERROR! Number of channels > max # of channels (%i/%i)\n",
            dec_noch, maxchannels);
    exit(0);
  }
  // END Header  **********************

  // Data Words  **********************
  for(int n=0; n< (int) dec_noch; n++){
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
    }
    dec_tow=((dataword>>24) & 0x7);  // decode type of dataword
    dec_data = (dataword & 0xfff);  // decode actual data
    if(maxchannels==32){ // ch on board max 32
      dec_ch =((dataword>>16) & 0x3f);  // decode channel number
    }
     else{ //ch on board max 16, N version of board
      dec_ch =((dataword>>17) & 0x1f);  // decode dataword
    }
    if(dec_ch>=maxchannels){
      printf("Decoded channel number is too big!\n");
      return -1;
    }
    if(dec_tow==0 && (dec_geo==(dataword>>27))){ // tow and geo is ok
      if(what==30 || what==31){ // tdc has a valid data bit
        if(((dataword>>14)&0x1)==0) dec_data=0; // valid data?
      }
      data[dec_ch] = dec_data;
    }
     else  printf("D");
  }
   // END Data Words  **********************

  // Trail  **********************
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  dec_tow=((dataword>>24) & 0x7); // decode type of word
  if(dec_tow!=4){
    printf("T"); //if it wasn't trail
    exit(0);
  }
  noe_of_board =((dataword) & 0xffffff);    //
   // END Trail **********************
// printf("END Reading out normal module %i\n", module_number);
  if(noe_of_board!=no_of_ev){  // if noe of board != noe external loop
    if(realign_first_noe==0){
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("A CAEN Module (%i) has got a wrong number of event\n", module_number);
      printf("no of event from module: %i, from counter %i\n", noe_of_board, no_of_ev);
      printf("EXIT\n");
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      return(2);
    }
    else{
     realign_first_noe=0;
     printf("Realigning NOE (%i) to %i\n", 
             module_number, no_of_ev-noe_of_board);
     int n;
     n=no_of_ev-noe_of_board;
     if(n>0) n*=-1;
     return(n);
    }
  }
  return(1);
}  // ReadRaw_class::read_one_event
    
unsigned int ReadRaw_class::get_noe(void){
  return noe_of_board;
}

int ReadRaw_class::get_maxchannels(void){
  return maxchannels;
}

int ReadRaw_class::get_maxevents(void){
  return max_event_per_channel;
}

int ReadRaw_class_BLANK::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  return 1;
}

void ReadRaw_class_BLANK::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  realign_first_noe=ext_realign_first_noe;
  module_number=n; what=ext_what;  in=ext_in;
  maxchannels=1; data = new unsigned int[1];
  read_further_parameters();
}

unsigned int ReadRaw_class_BLANK::getval(int channel, int event) const{
  printf("Blank, no data!\n"); return 0;
}


void ReadRaw_class_SIS3820::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  printf("Init module %i\n", n);
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=1;
  maxchannels=5;
  data = new unsigned int[5];
  read_further_parameters();
}  // ReadRaw_class::init

int ReadRaw_class_SIS3820::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word

  memset(data, 0, sizeof(data[0]));    // set all to 0 
  // header
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
//    printf("header sis: 0x%x\n", dataword);
  if(dataword!=0xaaaa){
    printf("Not a SIS3820 header!\n");
    exit(0);
  }
  // END header

  // data
  for(int n=0;n<5; n++){ //  counter 1-4, direct data in
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
    }
    data[n]= dataword;
//    printf("dataword sis: 0x%x\n", dataword);
  }
  // END data

  // trail
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  if(dataword!=0x5555){
    printf("Not a SIS3820 trail! (0x%x)\n", dataword);
    exit(0);
  }
  // END trail
  noe_of_board = data[1]-1;  // starting with 0

  if(noe_of_board!=no_of_ev){  // if noe of board != noe external loop
    if(realign_first_noe==0){
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("A SIS3820 Module (%i) has got a wrong number of event\n", module_number);
      printf("no of event from module: %i, from counter %i\n", noe_of_board, no_of_ev);
      printf("EXIT\n");
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      return(2);
    }
    else{
      realign_first_noe=0;
      printf("Realigning NOE (%i) to %i\n", 
             module_number, no_of_ev-noe_of_board);
      int n;
      n=no_of_ev-noe_of_board;
      if(n>0) n*=-1;
      return(n);
    }
  }
  return(1);
} //ReadRaw_class_SIS3820::read_one_event



void ReadRaw_class_SIS3302::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  printf("Init module %i\n", n);
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=1;
  maxchannels=8;

  read_further_parameters();

  SADC_N_of_samples=further_parameters[2];
  max_event_per_channel=SADC_N_of_samples;
  printf("Init new data array with %i events/ch\n", max_event_per_channel);
  data = new unsigned int[SADC_N_of_samples/2*8+8+2]; // raw head/trail/ch
}  // ReadRaw_class::init

int ReadRaw_class_SIS3302::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word

//  printf("SIS3302 readout %i\n", events);  
/*  if(events==0){ // before, SADC_N_of_samples, decoded in 
                   // further_parameters[2], is not known.
    SADC_N_of_samples=further_parameters[2];
    max_event_per_channel=SADC_N_of_samples;
    printf("Init new data array\n");
    delete data;
    data = new unsigned int[SADC_N_of_samples/2*8+8+2]; // raw head/trai/ch
  }
*/
  memset(data, 0, (SADC_N_of_samples/2*8+8+2)*4);    // set all to 0 
  // header
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
//    printf("header sis: 0x%x\n", dataword);
  if(dataword!=0xfadcfadc){
    printf("Not a SIS3302 header! (0x%x)\n", dataword);
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    printf("Dataword after: Not a SIS3302 header! (0x%x)\n", dataword);

    exit(0);
  }
  // END header

  // data
  for(int n=0;n<8; n++){ //  8 channels
      // channel header, (ch<<28)+SADC_N_of_samples;
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
              events, no_of_int, ((double) events)*100/ ((double)no_of_int), 
              file_number); 
    }

    for(int m=0; m<SADC_N_of_samples; m+=2){
      if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
      events++; 
      if(events%wie_oft==0){
        printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), 
           file_number); 
      }
//
      data[m/2+SADC_N_of_samples*n/2]= dataword;
//    printf("dataword sis: 0x%x\n", dataword);
    }
  }
//  printf("End of data loop\n");
  // END data
  // trail
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  if(dataword!=0xfadceee){ 
    printf("Not a SIS3302 trail! (0x%x)\n", dataword);
    exit(0);
  }
  // END trail
  noe_of_board = 0; //no_of_int;

  return(1);
} //ReadRaw_class_SIS3302::read_one_event

unsigned int ReadRaw_class_SIS3302::getval(int channel, int event) const{
  int pos;
  if(channel>=0 && channel<8 && event>=0 && event<SADC_N_of_samples){
    pos=(event/2)+(channel*SADC_N_of_samples/2);
    if(event%2==0) return((data[pos])&0xffff);
    else return((data[pos])>>16);
  }  
   else{
    printf("Channels out of range!\n");
    return 0;
  }
} //  ReadRaw_class_SIS3302::getval

unsigned int ReadRaw_class_SIS3302::get_trace(
                                   int channel, unsigned int *array) const{
  int pos;
  if(channel>=0 && channel<8){
    for(int n=0; n<SADC_N_of_samples; n++){
      pos=(n/2)+(channel*SADC_N_of_samples/2);
      if(n%2==0) array[n]=((data[pos])&0xffff);
       else array[n]=((data[pos])>>16);
    }
    return 1;
  }     
   else{
    printf("Channels out of range!\n");
    return 0;
  }
} //  ReadRaw_class_SIS3302::get_trace

int ReadRaw_class_SIS3302::get_raw_ch_enabled(int channel){
   return 1;
} //  ReadRaw_class_SIS3302::get_raw_ch_enabled

void ReadRaw_class_BONNSYNC::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  printf("Init module %i\n", n);
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=1;
  maxchannels=1;
  data = new unsigned int[1];
  read_further_parameters();
}  // ReadRaw_class::init

int ReadRaw_class_BONNSYNC::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word

  memset(data, 0, sizeof(data[0]));    // set all to 0 
  // header
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
//    printf("header bonnsync: 0x%x\n", dataword);
  if(dataword!=0xbbaa){
    printf("Not a BONNSYNC header!\n");
    exit(0);
  }
  // END header

  // data
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
         events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  data[1]= dataword;
//    printf("dataword bonnsync: 0x%x\n", dataword);

  // END data

 // trail
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  if(dataword!=0xaabb){
    printf("Not a BONNSYNC trail! (0x%x)\n", dataword);
    exit(0);
  }
  // END trail
  noe_of_board = data[1];  // starting with 0

  if(noe_of_board!=no_of_ev){  // if noe of board != noe external loop
    if(realign_first_noe==0){
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("A BONNSYNC Module (%i) has got a wrong number of event\n", module_number);
      printf("no of event from module: %i, from counter %i\n", noe_of_board, no_of_ev);
      printf("EXIT\n");
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      return(2);
    }
    else{
      realign_first_noe=0;
      printf("Realigning NOE (%i) to %i\n", 
             module_number, no_of_ev-noe_of_board);
      int n;
      n=no_of_ev-noe_of_board;
      if(n>0) n*=-1;
      return(n);
    }
  }

  return(1);
} //ReadRaw_class_BONNSYNC::read_one_event

void ReadRaw_class_BONNTRACK::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  printf("Init module %i\n", n);
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=10;
  maxchannels=384*6;
  data = new unsigned int[1];
  
  read_further_parameters();
}  // ReadRaw_class::init

int ReadRaw_class_BONNTRACK::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word
  unsigned int trig_count_silicon;
  int no_fiberhits, no_siliconhits;
  int geo;
  int height, channel, side, box;

  memset(data, 0, sizeof(data[0]));    // set all to 0 
  for(int n=0; n<(384); n++){ 
    siliconhits[n][0]=siliconhits[n][1]=siliconhits[n][2]=siliconhits[n][3]=0; 
  }
  for(int n=0; n<(16); n++){ fiberhits[n][0]=fiberhits[n][1]=0; }

  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  noe_of_board=dataword;  // 0


  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  trig_count_silicon=dataword; //1, nod 2


  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  no_fiberhits=dataword; //2, nod 3

  for(int n=0; n<no_fiberhits; n++){
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
    }      
    if(((dataword&0xF8000000) >> 27)==8){ // n+3  header with geo
      geo=(dataword&0x1f);   // geo, 0 or 1
      if(geo>1){
        printf("Wrong geo id(0x%x)!\n", dataword);
        exit(0);
      }
    }
    if(((dataword&0xF8000000) >> 27)==0){ // n+3 data
      channel=((dataword&0x3e00000)>>21);
      if(geo=0){
        fiberdata1[channel][ fiberhits[channel][geo] ]=dataword&0x1fffff;
        fiberhits[channel][geo]++;
      }
      if(geo=1){
        fiberdata2[channel][ fiberhits[channel][geo] ]=dataword&0x1fffff;
        fiberhits[channel][geo]++;
      }
    }
  }


  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
         events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }      
  no_siliconhits=dataword; // 3+fiberhits, nod 4+fiberhits

  for(int n=0; n<no_siliconhits; n++){
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
         events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
    }
                        // n+4+fiberhits,  nod n+5+fiberhits
    height  = ((dataword&0x00000fff));
    channel = ((dataword&0x001ff000) >>12);  // hopefully starting with 0
    side    = ((dataword&0x00200000) >>21);  //hopefully 1 or 2
    box     = ((dataword&0x03000000) >>24); // 1 or 2
    box++; // nope 0 and 1
    side++; // was not 1 or 2, but 0 or 1
    if(side==1 && box==1){
      silicondata_b1_s1[channel][siliconhits[channel][0]]=height;
      siliconhits[channel][0]++;
    }
    if(side==1 && box==2){
      silicondata_b1_s2[channel][siliconhits[channel][1]]=height;
      siliconhits[channel][1]++;
    }
    if(side==2 && box==1){
      silicondata_b2_s1[channel][siliconhits[channel][2]]=height;
      siliconhits[channel][2]++;
    }
    if(side==2 && box==2){
      silicondata_b2_s2[channel][siliconhits[channel][3]]=height;
      siliconhits[channel][3]++;
    }
  }
  if(noe_of_board!=no_of_ev){  // if noe of board != noe external loop
    if(realign_first_noe==0){
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("A BONNTRACK Module (%i) has got a wrong number of event\n", module_number);
      printf("no of event from module: %i, from counter %i\n", noe_of_board, no_of_ev);
      printf("EXIT\n");
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      return(2);
    }
    else{
      realign_first_noe=0;
      printf("Realigning NOE (%i) to %i\n", 
             module_number, no_of_ev-noe_of_board);
      int n;
      n=no_of_ev-noe_of_board;
      if(n>0) n*=-1;
      return(n);
    }
  }

  return(1);
} //ReadRaw_class_BONNTRACK::read_one_event

unsigned int ReadRaw_class_BONNTRACK::getval(int channel, int event) const{
  if(event<0 || event>=10) {
    printf("Events out of range!\n");
    return 0;
  }
  if(channel>=0 && channel < 16){ // fiber 1
    return(fiberdata1[channel][event]);
  }
  else if(channel>=16 && channel < 32){ // fiber 2
    return(fiberdata2[(channel-16)][event]);
  }
  else if(channel>=32 && channel < 416){ // silicon b1s1
    return(silicondata_b1_s1[(channel-32)][event]);
  }
  else if(channel>=416 && channel < 800){ // silicon b1s2
    return(silicondata_b1_s2[(channel-416)][event]);
  }
  else if(channel>=800 && channel < 1184){ // silicon b2s1
    return(silicondata_b2_s1[(channel-800)][event]);
  }
  else if(channel>=1184 && channel < 1568){ // silicon b2s2
    return(silicondata_b2_s2[(channel-1184)][event]);
  }
  else{ 
    printf("Channels out of range!\n");
    return 0;
  }
} //ReadRaw_class_BONNTRACK::getval


int ReadRaw_class_BONNTRACK::get_multihitno(int channel){
   if(channel>=0 && channel < 16){ // fiber 1
    return(fiberhits[channel][0]);
  }
  else if(channel>=16 && channel < 32){ // fiber 2
    return(fiberhits[channel][1]);
  }
  else if(channel>=32 && channel < 416){ // silicon b1s1
    return(siliconhits[channel][0]);
  }
  else if(channel>=416 && channel < 800){ // silicon b1s2
    return(siliconhits[channel][1]);
  }
  else if(channel>=800 && channel < 1184){ // silicon b2s1
    return(siliconhits[channel][2]);
  }
  else if(channel>=1184 && channel < 1568){ // silicon b2s2
    return(siliconhits[channel][3]);
  }
  else{ 
    printf("Channels out of range!\n");
    return 0; 
  }
} //  ReadRaw::get_multihit_no

void ReadRaw_class_AVM16::init(int n, FILE * ext_in, int ext_what,
                         int ext_realign_first_noe){
  printf("Init module %i\n", n);
  realign_first_noe=ext_realign_first_noe;
  module_number=n;
  what=ext_what;
  in=ext_in;
  max_event_per_channel=1;
  maxchannels=16;

  read_further_parameters();

  samples_to_readout=((further_parameters[10]+1)*2);

  for(int n=0; n<16; n++){  // which channels were read out raw?
    if(((further_parameters[4]>>n)&0x1)==0x1){  // if particular channel is enabled
      raw_ch_ena[n]=1; // this one is ena (for his)
//      if (dbg>2) printf("Enable rawdata: Bit %i of 0x%x is used\n",
      //                       n, cha_raw);
    }
    else raw_ch_ena[n]=0;
  }

  max_event_per_channel=samples_to_readout;
  printf("Init new data array with %i events/ch\n", max_event_per_channel);
  data = new unsigned int[68*1024+2*4]; 
                                      // 68k max buffer + raw head/trail/ch
}  // ReadRaw_class_AVM16::init

int ReadRaw_class_AVM16::read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev){
  unsigned int dataword; // current data word
  int noe;
  int card_address, channel, data_ident;

  memset(data, 0, (68*1024+2*4)*4);    // set all to 0 
  memset(integral, 0, 16*6*4);
  memset(time, 0, 16*6*4);
  memset(events_integral, 0, 16*4);
  memset(events_time, 0, 16*4);
  memset(rawdatanumber, 0, 16*4);

  // header
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  if(dataword!=0xfadcfafa){
    printf("Not a AVM16 header! (0x%x)\n", dataword);
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    printf("Dataword after: Not a AVM16 header! (0x%x)\n", dataword);

    exit(0);
  }
  // END header

  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  eventlen=(dataword & 0x7fffffff)/4;

  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  // ignore timestamp

  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  noe_of_board = dataword; 
  
  // fill data[], integral[16][6], time[16][6], events_time[16], events_integral[16];

  for(int n=0;n<(eventlen-3); n++){ //  8 channels
    if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
    events++; 
    if(events%wie_oft==0){
      printf("%i of %i read in (%3.0f%%) from file number %i\n",
         events, no_of_int, ((double) events)*100/ ((double)no_of_int), 
         file_number); 
    }
    card_address = ((dataword>>28)&0xf);
    channel = ((dataword>>22)&0xf);
    data_ident = ((dataword>>16)&0x3f);

    if(((dataword>>12)&0x3ff)==0x0){  // raw data
      if(rawdatanumber[channel]<samples_to_readout){
        // copy rawdata to his
        data[channel*samples_to_readout+rawdatanumber[channel]]= (int) (dataword&0xfff);
      }
      rawdatanumber[channel]++;
    }
    else{
      if(((data_ident>>4)&0x3)==0x2){ // Integral
	if(events_integral[channel]<6){ // only first 6 integrals are allowed
  	  integral[channel][events_integral[channel]]= (int) dataword&0xffff;  // really 20bit?
	  events_integral[channel]++;
	}
      }
      else switch(data_ident){
	case 0x30: // window start time
	  if(events_time[channel]<6){ // only first 6 integrals are allowed
            time[channel][events_time[channel]]= (int) dataword&0xffff;
	    events_time[channel]++;
	  }
          break;
        case 0x31:
        case 0x32:
        case 0x33:
        case 0x34:
        case 0x35:
	  case 0x3d: // overflow bit
        case 0x36:
        case 0x37:
        case 0x20:
        case 0x43:
	  break;
        default:
	  printf("Error: Unknown data_ident: 0x%x (module %i)\n",
                  data_ident, module_number);
	  exit(0);
//    printf("dataword sis: 0x%x\n", dataword);
      }
    }
  }
//  printf("End of data loop\n");
  // END data

  // trail
  if(0==fread((char*) &dataword, sizeof(dataword), 1, in)) return(0); 
  events++; 
  if(events%wie_oft==0){
    printf("%i of %i read in (%3.0f%%) from file number %i\n",
           events, no_of_int, ((double) events)*100/ ((double)no_of_int), file_number); 
  }
  if(dataword!=0xcdafcdaf){ 
    printf("Not a AVM16 trail! (0x%x)\n", dataword);
    exit(0);
  }
  // END trail

  if(noe_of_board!=no_of_ev){  // if noe of board != noe external loop
    if(realign_first_noe==0){
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("A AVM16 Module (%i) has got a wrong number of event\n", module_number);
      printf("no of event from module: %i, from counter %i\n", noe_of_board, no_of_ev);
      printf("EXIT\n");
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
      return(2);
    }
    else{
      realign_first_noe=0;
      printf("Realigning NOE (%i) to %i\n", 
             module_number, no_of_ev-noe_of_board);
      int n;
      n=no_of_ev-noe_of_board;
      if(n>0) n*=-1;
      return(n);  
    }
  }
  return(1);
} //ReadRaw_class_AVM16::read_one_event

unsigned int ReadRaw_class_AVM16::getval(int channel, int event) const{
  int pos;
  if(channel>=0 && channel<16 && event>=0 && event<samples_to_readout){
    pos=(event)+(channel*samples_to_readout);
    return((data[pos])>>16);
  }  
  else if(channel>=16 && channel<32 && event>=0 && event<6){ // feature extr. integral
    return(integral[channel-16][event]);
  }
  else if(channel>=32 && channel<48 && event>=0 && event<6){ // feature extr. time
    return(time[channel-32][event]);
  }
  else{
    printf("Channels out of range!\n");
    return 0;
  }
} //  ReadRaw_class_AVM16::getval

unsigned int ReadRaw_class_AVM16::get_trace(
                                   int channel, unsigned int *array) const{
  int pos;
  if(channel>=0 && channel<16){
    for(int n=0; n<samples_to_readout; n++){
      pos=(n)+(channel*samples_to_readout);
      array[n]=(data[pos]);
    }
    return 1;
  } 
   else{
    printf("Channels out of range!\n");
    return 0;
  }
} //  ReadRaw_class_AVM16::get_trace

int ReadRaw_class_AVM16::get_multihitno(int channel){
  if(channel>=0 && channel < 16){ // raw data
    return(1);
  }
  else if(channel>=16 && channel < 32){ // integral
    return(events_integral[channel]);
  }
  else if(channel>=32 && channel < 48){ // time
    return(events_time[channel]);
  }
  else if(channel>=48 && channel < 64){ // own_ph
    return(1);
  }
  else if(channel>=64 && channel < 80){ // pedestal
    return(1);
  }
  else{
    printf("Unknown channel %i\n", channel);
    exit(0);
  }
}

int ReadRaw_class_AVM16::get_raw_ch_enabled(int channel){
  if(channel<16){  
    return raw_ch_ena[channel];
  }
  else{
    printf("AVM16 channel out of bound\n");
    exit(0);
  }
} //  ReadRaw_class_AVM16::get_raw_ch_enabled
