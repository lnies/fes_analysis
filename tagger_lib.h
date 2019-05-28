#ifndef _LIB_TAGGER_H
#define _LIB_TAGGER_H

struct dataheaderstruct{
	unsigned int a;
	unsigned int b;
	unsigned int EventCB;
	unsigned int EventProto;
	unsigned int nTaggerhits;
	unsigned int nCherenkovhits;
}; 

class Taggerfile_class{
  protected:
    static int initialized;
    FILE * in;
  public:
    Taggerfile_class(){};
    ~Taggerfile_class(){};
    void init(char * file);
    int read_one_event(unsigned int no_of_ev);
    double get_taggertime(int number);
    double get_taggerenergy(int number);
    double get_cherenkovtime(int number);
    int get_no_taggerhits(void);
    int get_no_cherenkovhits(void);

    dataheaderstruct dataheader;
    double beamphotons[32][2];  // 0 energy, 1 time
    double cherenkovhits[32];
};

#endif
