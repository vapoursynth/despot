/*
    DeSpot - Conditional Temporal Despotting Filter for Avisynth 2.5

    Include file

    Version 3.5  November 17, 2005.

Copyright (C)2003-2005 Alexander G. Balakhnin aka Fizick.
bag@hotmail.ru         http://bag.hotmail.ru
under the GNU General Public Licence version 2.

This plugin is based on Conditional Temporal Median Filter (c-plugin) for Avisynth 2.5
Version 0.93, September 27, 2003
Copyright (C) 2003 Kevin Atkinson (kevin.med@atkinson.dhs.org) under the GNU GPL version 2.
http://kevin.atkinson.dhs.org/temporal_median/

*/
#include <string>
#include <avs/types.h>

#define M_SEARCH 0
#define M_MEDIAN 4

#define S_DENOISE 0
#define S_MARK    1
#define S_MAP     2

struct Parms {
    size_t size;
    int width, height;
    int pitch;
    BYTE p1 = 24;
    BYTE p2 = 12;
    short pwidth = 16;
    short pheight = 5;
    int y_next = 1;
    BYTE mthres = 16;
    int mwidth = 7;
    int mheight = 5;
    int merode = 33;
    bool median = false;
    int show = 0;
    BYTE mark_v = 255;
    bool show_chroma = false;
    bool ranked = true;
    int sign = 0;
    int maxpts = 0;
    int p1percent = 10;
    int dilate = 1;
    bool fitluma = true;
    BYTE blur = 1;
    int tsmooth = 0;
    bool motpn = true;
    int seg = 2;
    bool color = false;
    int mscene = 40;
    int minpts = 0;

    std::string	outfilename;
    bool mc_flag;
    int spotmax1;
    int spotmax2;
};

constexpr char B_NOTHING = 0;
constexpr char B_SEGMENT = 1; // single segment
constexpr char B_LINK = 2; // has master
constexpr char B_MASTER = 4; // has slave link

struct Segment {
    bool p1yes = false;
    char what = B_SEGMENT;
    short ly; // top
    short lx1; // left
    short lx2; // right
    long nump1 = 0;
    long nump2 = 0;
    long sumnoise = 0;
    struct Data {
        struct BStruct {
            short x1;
            short x2;
            short y1;
            short y2;
        };
        struct SStruct {
            short width;
            short height;
            bool is_noise;
        };
        BStruct b = {};
        SStruct s = {};
        Segment *nowis = nullptr;

        explicit Data(short x, short y) {
            b.x1 = x;
            b.x2 = x;
            b.y1 = y;
            b.y2 = y;
        }
    } data;

    Segment(short ly0, short lx10) :
        ly(ly0), lx1(lx10), lx2(lx10), data(lx10, ly0) {
    }
};

struct Segments {
    Segment *data;
    Segment *end;
    Segments() : data(0), end(0) {
    }
    void clear() {
        end = data;
    }
};

typedef void Exec(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch,
    BYTE *c_noise, BYTE *c_motion,
    const BYTE *n, int Npitch,
    BYTE *o, int Opitch,
    const Parms &);// BYTE * n_motion, - remove parameter in v.3.0



void find_sizes(BYTE *noise, Segments &segs, const Parms &); // change parameters order in v.3.0

void mark_noise(Segments &segs, BYTE *noise, const Parms &);

void find_motion(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    BYTE *motion, const Parms &);
Exec cond_median;

void motion_denoise(BYTE *moving, const Parms &); // 3.22

void find_outliers(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    const BYTE *n, int Npitch, BYTE *c_noise, const Parms &);

void mark_motion(const BYTE *p, int Ppitch, BYTE *p_noise,
    const BYTE *c, int Cpitch, BYTE *c_noise, BYTE *c_motion,
    const Parms &);

Exec remove_outliers;
Exec mark_outliers;
Exec map_outliers;

void noise_dilate(BYTE *noise, const Parms &); // added in v.1.3
void reject_on_motion(Segments &segs, BYTE *motion, const Parms &p); // added in v.3.0
void motion_merge(BYTE *c_motion, BYTE *n_motion, BYTE *m_motion, const Parms &parms); // added in v.3.0
void noise_to_one(BYTE *p_noise, BYTE *c_noise, BYTE *c_motion, const Parms &parms);// added in v.3.0
void remove_segments(Segments &segs, const BYTE *P, int Ppitch, const BYTE *C, int Cpitch, const BYTE *N, int Npitch, BYTE *o, int Opitch, BYTE *c_noise, const Parms &parms);// added in v.3.0
void temporal_smooth(Segments &segs, const BYTE *P, int Ppitch, const BYTE *N, int Npitch, BYTE *o, int Opitch, BYTE *c_noise, BYTE *c_motion, const Parms &parms);// added in v.3.0 , changed in 3.2
void clean_color_plane(const BYTE *p, int ppitch, const BYTE *c, int cpitch, const BYTE *n, int npitch, int width, int height, BYTE *c_noise, BYTE *o, int opitch, Parms &Params);// added in v.3.1
void mark_color_plane(BYTE *c_noise, BYTE *oV, int opitchV, int widthUV, int heightUV, Parms &Params);// added in v.3.1
void motion_scene(BYTE *motion, const Parms &parms); // added in v.3.3
void add_external_mask(const BYTE *ext, int ext_pitch, BYTE *motion, int pitch, int width, int height); // added in v.3.5
void print_segments(const Segments &segs, FILE *f_ptr, int fn, double frate, int w, int h, int spotmax1, int spotmax2);
