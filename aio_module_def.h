#ifndef _AIO_MODULE_DEF_H
#define _AIO_MODULE_DEF_H

#define NTEC_MULTI 11
#define VETO_MULTI 21
#define V775N_TDC 30
#define V775_TDC 31
#define V785N_PSADC 40
#define V785_PSADC  41
#define V792N_QDC  50
#define V792_QDC   51
#define V965N_DRQDC 60
#define V965_DRQDC  61
#define IOL  89
#define SIS3820_SYNC 90
#define BONN_SYNC 91
#define V895_LED 110
#define AMV16_FLADC 150
#define SIS3302_FADC 160
#define BONN_TRACK 170
#define AVM16_FADC 180

#define AIO_max_no_of_boards 16
#define AIO_versionsstr "v11.0.0"
#define AIO_default_filename "data_all_v11"


void AIO_what_to_name(int what, char *name);
void timestring(char * string);
void binarystring(int value, char* string);


#endif
