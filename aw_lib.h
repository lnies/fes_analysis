#ifndef _LIB_AWALLINONE_H
#define _LIB_AWALLINONE_H

#include "aio_module_def.h"


class ReadRaw_class{
  protected:
    static unsigned int wie_oft;
    static unsigned int events;
    FILE * in;
    int max_no_rawdata;
    int adc_range;
    int module_number;
    int what;
    int maxchannels;
    int max_event_per_channel;
    int further_parameters[20];  // in  init_module_variables
    int realign_first_noe;
    unsigned int noe_of_board;
    unsigned int *data;
  public:
    ReadRaw_class(){
    }
    virtual ~ReadRaw_class(){ 
      printf("Destructor!\n");
      delete [] data;
      fclose(in);
    }
    virtual void init(int n, FILE * ext_in, int ext_what, 
                      int ext_realign_first_noe);
    void set_filepointer(FILE * ext_in);
    virtual unsigned int getval(int channel, int event) const;
    virtual int get_multihitno(int channel);
    virtual unsigned int get_trace(int channel, unsigned int *array) const;
    virtual int read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev);
    virtual int get_raw_ch_enabled(int channel);
    void read_further_parameters(void);
    unsigned int get_noe(void);
    int get_maxchannels(void);
    int get_maxevents(void);
    int get_further_parameter(int number);
};

class ReadRaw_class_BLANK: public ReadRaw_class{
  public:
    virtual int read_one_event(unsigned int no_of_int, int file_number, unsigned int no_of_ev);
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual unsigned int getval(int channel, int event) const;
};

class ReadRaw_class_SIS3820: public ReadRaw_class{
  public:
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual int read_one_event(unsigned int no_of_int, int filenumber, unsigned int no_of_ev);
};

class ReadRaw_class_SIS3302: public ReadRaw_class{
  protected:
    int SADC_N_of_samples;
  public:
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual int read_one_event(unsigned int no_of_int, int filenumber, unsigned int no_of_ev);
    virtual unsigned int getval(int channel, int event) const;
    virtual unsigned int get_trace(int channel, unsigned int *array) const;
    virtual int get_raw_ch_enabled(int channel);
};

class ReadRaw_class_BONNSYNC: public ReadRaw_class{
  public:
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual int read_one_event(unsigned int no_of_int, int filenumber, unsigned int no_of_ev);
};

class ReadRaw_class_BONNTRACK: public ReadRaw_class{
 protected:
    int siliconhits[384][4];
    int fiberhits[16][2];
    int fiberdata1[16][10];
    int fiberdata2[16][10];
    int silicondata_b1_s1[384][10];
    int silicondata_b1_s2[384][10];
    int silicondata_b2_s1[384][10];
    int silicondata_b2_s2[384][10];
  public:
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual int read_one_event(unsigned int no_of_int, int filenumber, unsigned int no_of_ev);
    virtual unsigned int getval(int channel, int event) const;
    virtual int get_multihitno(int channel);
};

class ReadRaw_class_AVM16: public ReadRaw_class{
  protected:
    int integral[16][6];
    int time[16][6];
    int events_time[16], events_integral[16];
    int raw_ch_ena[16];
    int samples_to_readout;
    int rawdatanumber[16];
    int eventlen;
  public:
    virtual void init(int n, FILE * ext_in, int ext_what,
                      int ext_realign_first_noe);
    virtual int read_one_event(unsigned int no_of_int, int filenumber, unsigned int no_of_ev);
    virtual unsigned int getval(int channel, int event) const;
    virtual unsigned int get_trace(int channel, unsigned int *array) const;
    virtual int get_multihitno(int channel);
    virtual int get_raw_ch_enabled(int channel);
};



class ReadSystem_class{
  protected:
    char filename[100];
    char file[100];
    int single_file;
    int file_number;
    FILE * in;
    int number_of_boards;
    int max_no_rawdata[AIO_max_no_of_boards];
    int no_of_boards, what[16], geos[16];
    int old_no_of_boards, old_what[16], old_geos[16];
    int no_of_int;
    char gen_starttime[50];
    char starttime[50];
    char stoptime[50];
    unsigned int geo_to_bnr[50];
    int realign_first_noe;
    ReadRaw_class *read_modules[16];
  public:
    ReadSystem_class(){
    }
    ~ReadSystem_class(){
    }
    int init_readout(char * file, int ext_realign_first_noe);
    int read_one_event(unsigned int & no_of_ev);
    unsigned int get_value(int board, int channel, int event);
    unsigned int get_trace(int board, int channel, unsigned int *array);
    unsigned int get_noe(int module);
    int get_maxchannels(int module);
    int get_maxevents(int module);
    int get_multihitno(int module, int channel);
    int get_further_parameter(int module, int number);
    int get_raw_ch_enabled(int module, int channel);
  private:
    int open_file(void);	
    void read_in_header(void);
};

#endif
