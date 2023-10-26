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
#include <VapourSynth4.h>
#include <VSHelper4.h>

constexpr int S_DENOISE = 0;
constexpr int S_MARK = 1;
constexpr int S_MAP = 2;

struct Parms {
    size_t size = 0;
    int width = -1;
    int height = -1;
    int pitch = -1;
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
    bool mc_flag = false;
    int spotmax1 = 12;
    int spotmax2 = 20;
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

        Data(short x, short y) {
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
    Segment *data = nullptr;
    Segment *end = nullptr;
    void clear() {
        end = data;
    }
};

typedef void Exec(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch,
    BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion,
    const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch,
    const Parms &);

void find_sizes(BYTE *VS_RESTRICT noise, Segments &segs, const Parms &);

void mark_noise(Segments &segs, BYTE *VS_RESTRICT noise, const Parms &);

void find_motion(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    BYTE *VS_RESTRICT motion, const Parms &);
Exec cond_median;

void motion_denoise(BYTE *VS_RESTRICT moving, const Parms &);

void find_outliers(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    const BYTE *n, int Npitch, BYTE *VS_RESTRICT c_noise, const Parms &);

Exec remove_outliers;
Exec mark_outliers;
Exec map_outliers;

void noise_dilate(BYTE *VS_RESTRICT noise, const Parms &);
void reject_on_motion(Segments &segs, BYTE *VS_RESTRICT motion, const Parms &p);
void motion_merge(BYTE *VS_RESTRICT c_motion, BYTE *VS_RESTRICT n_motion, BYTE *VS_RESTRICT m_motion, const Parms &parms);
void noise_to_one(BYTE *VS_RESTRICT p_noise, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion, const Parms &parms);
void remove_segments(Segments &segs, const BYTE *P, int Ppitch, const BYTE *C, int Cpitch, const BYTE *N, int Npitch, BYTE *VS_RESTRICT o, int Opitch, BYTE *VS_RESTRICT c_noise, const Parms &parms);
void temporal_smooth(Segments &segs, const BYTE *P, int Ppitch, const BYTE *N, int Npitch, BYTE *VS_RESTRICT o, int Opitch, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion, const Parms &parms);
void clean_color_plane(const BYTE *p, int ppitch, const BYTE *c, int cpitch, const BYTE *n, int npitch, int width, int height, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT o, int opitch, Parms &Params);
void mark_color_plane(BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT oV, int opitchV, int widthUV, int heightUV, Parms &Params);
void motion_scene(BYTE *VS_RESTRICT motion, const Parms &parms);
void add_external_mask(const BYTE *ext, int ext_pitch, BYTE *VS_RESTRICT motion, int pitch, int width, int height);
void print_segments(const Segments &segs, FILE *f_ptr, int fn, double frate, int w, int h, int spotmax1, int spotmax2);
