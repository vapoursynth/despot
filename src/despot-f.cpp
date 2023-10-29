/*
    DeSpot - Conditional Temporal Despotting Filter for Avisynth 2.5

    Filter functions file
    Version 3.5 November 25, 2005.


Copyright (C)2003-2004 Alexander G. Balakhnin aka Fizick.
bag@hotmail.ru         http://bag.hotmail.ru
under the GNU General Public Licence version 2.

This plugin is based on Conditional Temporal Median Filter (c-plugin) for Avisynth 2.5
Version 0.93, September 27, 2003
Copyright (C) 2003 Kevin Atkinson (kevin.med@atkinson.dhs.org) under the GNU GPL version 2.
http://kevin.atkinson.dhs.org/temporal_median/

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "despot.hpp"

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <malloc.h>

static void memzero(void *d0, size_t s) {
    memset(d0, 0, s);
}

static void combine_xs(BYTE *VS_RESTRICT noise, Segments &segs, const Parms &p) {
    BYTE *l = noise;
    Segment *c = nullptr;

    BYTE p2 = p.p2;

    for (int y = 0; y != p.height; ++y) {
        for (int x = 0; x != p.width; ++x) {
            if (l[x] > p2) {

                if (c == 0) { // new left edge found
                    c = segs.end++; // get next free segment number
                    *c = Segment(y, x); // fixme, why was this so roundabout with memcpy?
                }

                c->nump2 += 1;
                c->sumnoise += l[x];

                if (l[x] > p.p1) {
                    c->p1yes = true;
                    c->nump1 += 1;
                }

            } else if (c != 0) { // right edge found

                c->data.b.x2 = x - 1;
                c->lx2 = x - 1;
                c = 0;

            }
        }
        if (c != 0) {
            c->lx2 = p.width - 1;
            c = 0;
        }
        l += p.pitch;
    }
}


void combine_ys(Segments &segs, const Parms &p) {
    segs.end->ly = p.height; // to avoid some special cases
    Segment *c = segs.data; // first segment
    Segment *n = c + 1;
    while (n->ly == (n - 1)->ly) ++n; // search next segment with ly not equal to ly of c
    Segment *n0;

    bool c_has_next, n_has_next;

    while (n < segs.end) {
        n0 = n;
        if (n->ly < c->ly + p.y_next) // search next segment with ly > ly of c
        {
            ++n; while (n->ly == (n - 1)->ly) ++n;
        }

        if (c->ly + p.y_next == n->ly) { // we found next segment in neighbor line
            do {
                if ((c->lx1 <= n->lx1 && c->lx2 >= n->lx1 - 1) ||
                    (c->lx1 > n->lx1 && n->lx2 >= c->lx1 - 1)) // overlap or neighbor
                {
                    Segment *a = c; while (a->what & B_LINK) a = a->data.nowis;
                    Segment *b = n; while (b->what & B_LINK) b = b->data.nowis;
                    if (a > b) {
                        Segment *t = a; a = b; b = t;
                    } // this is important
                    a->p1yes = a->p1yes || b->p1yes; // Changed par name in v 1.2
                    a->data.b.x1 = std::min(a->data.b.x1, b->data.b.x1);
                    a->data.b.x2 = std::max(a->data.b.x2, b->data.b.x2);
                    a->data.b.y2 = n->ly; // and what about y1 ? - not changed


                    if (a != b) {
                        a->nump1 += b->nump1;
                        a->nump2 += b->nump2;
                        a->sumnoise += b->sumnoise;
                        a->what = a->what | B_MASTER;
                        b->what = b->what | B_LINK;
                        b->data.nowis = a; // "b" now is included in "a" as LINK
                    }
                }

                c_has_next = c->ly == (c + 1)->ly;
                n_has_next = n->ly == (n + 1)->ly;
                if (c->lx2 <= n->lx2 && c_has_next) ++c;
                else if (n_has_next)                     ++n;
                else if (c_has_next)                     ++c;
                else break;

            } while (true);

            ++c;
            ++n;

        } else {
            if (c->ly < n0->ly && n0->ly < n->ly) {
                c = n0;
            } else {
                c = n;
                ++n;
                while (n->ly == (n - 1)->ly) ++n;
            }
        }
    }
}

static void to_size(Segments &segs, const Parms &p) {

    const short PH = p.pheight * p.y_next;
    for (Segment *cur = segs.data; cur != segs.end; ++cur) {
        if (!(cur->what & B_LINK)) {  // use only super-master or single segment - rewrited in v.3.0
            Segment::Data &b = cur->data;
            short w = b.b.x2 - b.b.x1 + 1;
            short h = b.b.y2 - b.b.y1 + 1;
            b.s.is_noise = cur->p1yes && w <= p.pwidth && h <= PH;

            if (p.maxpts > 0 && b.s.is_noise) {
                // accept only  segments with not very many points
                b.s.is_noise = cur->nump2 <= p.maxpts;
            }

            if (p.minpts > 0 && b.s.is_noise) {
                // accept only  segments with not very many points
                b.s.is_noise = cur->nump2 >= p.minpts;
            }

            if (p.p1percent > 0 && b.s.is_noise) {
                // accept only strong segments with big relative number of truth spot (not outliers)
                b.s.is_noise = (100 * cur->nump1) / cur->nump2 >= p.p1percent;
            }

            b.s.width = w; // where are it used ?
            b.s.height = h; // And it will also change x1 x2 y1 y2 nowis in union? I changed union to structure
        } else { //  B_LINK
            // set data is_noise as equal of data is_noise of master segment (w, h, is_noise) - changed in v.3.0
            // Now we not lost info about master (nowis)
            Segment *master = cur->data.nowis;
            cur->data.s.is_noise = master->data.s.is_noise;
        }
    }
}

// scan noise and create segments

void find_sizes(BYTE *VS_RESTRICT noise, Segments &segs, const Parms &parms) {
    combine_xs(noise, segs, parms);
    combine_ys(segs, parms);
    to_size(segs, parms);
}

void mark_noise(Segments &segs, BYTE *VS_RESTRICT l, const Parms &p) {
    memzero(l, p.size);
    Segment *c = segs.data;
    segs.end->ly = p.height; // to avoid having to check for < segs.end
    for (int y = 0; y != p.height; ++y) {
        for (; c->ly == y; ++c) {
            if (c->data.s.is_noise) {
                BYTE *b = l + c->lx1;
                BYTE *end = l + c->lx2 + 1;
                for (; b != end; ++b) *b = 0xFF; // 255 is noise (spots)
            }
        }
        l += p.pitch;
    }
}

//
// Reject segments near motion
//


void reject_on_motion(Segments &segs, BYTE *VS_RESTRICT motion, const Parms &p) {
    Segment *master;
    Segment *mastertmp;

    // reverse order is more fast because most segments has master segment
    Segment *c = segs.end - 1;
    BYTE *mot = motion + (p.height - 1) * p.pitch;

    // First pass - search segments over or near motion pixels and reject its  and all its masters
    for (int y = p.height - 1; y >= 0; y--) {

        while (c >= segs.data) { // next segment with y
            if (c->ly == y) {
                bool is_noise = c->data.s.is_noise;
                if (is_noise) { // check current line y
                    BYTE *b = mot + std::max(c->lx1 - 1, 0); // left edge-1 to check neighbor points for motion
                    BYTE *end = mot + std::min(c->lx2 + 2, p.width); // right edge+1
                    for (; b < end; ++b) {
                        if (*b != 0) { // if found motion pixel on or near spot(segment)
                            c->data.s.is_noise = false; // reject this segment
                            if (p.seg == 2) { // 2D spot mode - mark master segments
                                master = c->data.nowis; // get master segment
                                while (master != 0) { // check if master exist
                                    master->data.s.is_noise = false; // reject master
                                    mastertmp = master; // copy
                                    master = mastertmp->data.nowis; // get higher master
                                }
                            }
                            break; // exit this for loop and go to next segment
                        }
                    }

                    if (y >= p.y_next) {// also check vertical motion neighbors at line y-1
                        BYTE *b = mot + std::max(c->lx1 - 1, 0) - p.pitch * p.y_next; // left edge-1 to check neighbor points for motion
                        BYTE *end = mot + std::min(c->lx2 + 2, p.width) - p.pitch * p.y_next; // right edge+1
                        for (; b < end; ++b) {
                            if (*b != 0) { // if found motion pixel on or near spot(segment)
                                c->data.s.is_noise = false; // reject this segment
                                if (p.seg == 2) {// 2D spot mode- mark master segments
                                    master = c->data.nowis; // get master segment
                                    while (master != 0) { // check if master exist
                                        master->data.s.is_noise = false; // reject master
                                        mastertmp = master; // copy
                                        master = mastertmp->data.nowis; // get higher master
                                    }
                                }
                                break; // exit this for loop and go to next segment
                            }
                        }
                    }

                    if (y < p.height - p.y_next) {// also check vertical motion neighbors at line y+1
                        BYTE *b = mot + std::max(c->lx1 - 1, 0) + p.pitch * p.y_next; // left edge-1 to check neighbor points for motion
                        BYTE *end = mot + std::min(c->lx2 + 2, p.width) + p.pitch * p.y_next; // right edge+1
                        for (; b < end; ++b) {
                            if (*b != 0) { // if found motion pixel on or near spot(segment)
                                c->data.s.is_noise = false; // reject this segment
                                if (p.seg == 2) {  // 2D spot mode- mark master segments
                                    master = c->data.nowis; // get master segment
                                    while (master != 0) { // check if master exist
                                        master->data.s.is_noise = false; // reject master
                                        mastertmp = master; // copy
                                        master = mastertmp->data.nowis; // get higher master
                                    }
                                }
                                break; // exit this for loop and go to next segment
                            }
                        }
                    }

                }
            } else 	break;

            c--; //next prev
        }
        mot -= p.pitch;
    }

    if (p.seg == 2) {  // 2D spot mode
        // second pass - reject all slaves of rejected masters
        for (c = segs.data; c != segs.end; ++c) {
            master = c->data.nowis;
            if (master != 0)
                c->data.s.is_noise = master->data.s.is_noise;
        }
    }

}


//
// Motion Detection
//

static void find_mot_line(const BYTE *p, const BYTE *c, BYTE *VS_RESTRICT motion, int width, int mthres) {
    for (int x = 0; x < width; ++x) {
        if (abs(p[x] - c[x]) > mthres)
            motion[x] = 3; //  motion without noise
        else
            motion[x] = 0; // 0 is no motion
    }
}


void find_motion(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch, BYTE *VS_RESTRICT motion,
    const Parms &parms) {
    int width = parms.pitch;
    int mthres = parms.mthres;
    for (int y = 0; y != parms.height; ++y) {
        find_mot_line(p, c, motion, width, mthres);

        p += Ppitch;
        c += Cpitch;
        motion += parms.pitch;
    }
}

static void motion_summation(const BYTE *moving, signed char *VS_RESTRICT r, int w, int h, int y1, int y2) {
    if (y1 >= 0) {
        const BYTE *r2 = moving + w * y1;
        for (int x = 0; x < w; ++x) {
            r[x] -= r2[x];
        }
    }
    if (y2 < h) {
        const BYTE *r2 = moving + w * y2;
        for (int x = 0; x < w; ++x) {
            r[x] += r2[x];
        }
    }
}

void motion_denoise(BYTE *VS_RESTRICT moving, const Parms &p) {
    // This function is based on the motion map denoise code found in
    // Donald Graft's Smart Deinterlacer Filter.  I (K.A.) completely rewrote
    // it to so that its complexity is O(w*h) instead of
    // O(w*h*MWIDTH*MHEIGHT).

      // comments  v.1.3:
      // I do not understand this code in detail.
      // Erode stage removes small objects of motion map,
      // dilate stage enlarge remaining objects.
      // seems,at input, moving[x] may be
      //=3 ( = mratio in v.2.1 with 3 default) if only motion, or
      //=1 if motion and noise exist at this pixel
      // at output, moving[x] may be 0 (no motion) or 255 (if motion)

    const int w = p.pitch;
    const int h = p.height;
    const int w_8 = w / 8;

    const int MHEIGHT_2 = p.mheight / 2;
    const int MWIDTH_2 = p.mwidth / 2;
    const int Y_NEXT = p.y_next;
    const int MPx3 = (p.merode * 3 * (p.mwidth + 1) * (p.mheight + 1)) / 100; // multiplied by mratio (default =3) because motion pixels mark=mratio


    if (MHEIGHT_2 == 0 && MWIDTH_2 == 0)
        return;

    //(only dynamic arrays with non-constant size):
    BYTE *fmoving = (BYTE *)malloc(sizeof(BYTE) * h * w);

    // (only dynamic arrays with non-constant size):
    size_t rowpadded_size = ((w + 64) / 8) * 8 + 8;
    signed char *row = (signed char *)malloc(sizeof(char) * p.y_next * rowpadded_size);

    if (MPx3 > 0) {
        // Erode

        for (int i = 0; i < Y_NEXT; i++) {
            memzero(&row[i * rowpadded_size], rowpadded_size);
        }

        for (int y = -(MHEIGHT_2)*Y_NEXT; y < h; ++y) {
            signed char *r = row + (y & (Y_NEXT - 1)) + 32;
            int y1 = y - (MHEIGHT_2)*Y_NEXT - Y_NEXT;
            int	y2 = y + (MHEIGHT_2)*Y_NEXT;
            motion_summation(moving, r, w, h, y1, y2);


            if (y >= 0) {
                BYTE *r2 = moving + w * y;
                BYTE *fw = fmoving + w * y;
                int total = 0;
                for (int x = 0; x != MWIDTH_2; ++x)
                    total += r[x]; //  summa of motion pixels marks (but not simple! remember about mratio(3) and 1)
                for (int x = 0; x != w; ++x) {
                    total -= r[x - MWIDTH_2 - 1];
                    total += r[x + MWIDTH_2];
                    fw[x] = (r2[x] && (total > MPx3)) ? 1 : 0;
                }
            }
        }
        // comments :
          // fmoving[x] is intermediate output after  erode
          // seems, fmoving[x]=1, if current pixel x have motion,
        // and total number of its neighbors at region MWidth+1 x MHeight+1 > MP
        // may be noise pixel also give some addition (not 3 but 1)
    } else { // no erode (faster)
        // we do not use erode, but we must initialize fmoving,
          //simply from  motion array data, but of 1 and 0 (not 3 and 0)
          // code added instead of erode
        for (int y = 0; y < h; y++) {
            BYTE *r2 = moving + w * y; // lines
            BYTE *fw = fmoving + w * y;
            for (int x = 0; x < w; x++) {
                fw[x] = r2[x] ? 1 : 0;
            }
        }

    }
    // Dilate
    for (int i = 0; i < Y_NEXT; i++) {
        memzero(&row[i * rowpadded_size], rowpadded_size);
    }

    for (int y = -(MHEIGHT_2)*Y_NEXT; y < h; ++y) {
        signed char *r = row + (y & (Y_NEXT - 1)) + 32;
        int y1 = y - (MHEIGHT_2)*Y_NEXT - Y_NEXT;
        int y2 = y + (MHEIGHT_2)*Y_NEXT;

        motion_summation(fmoving, r, w, h, y1, y2);

        if (y >= 0) {
            BYTE *fw = moving + w * y;
            int total = 0;
            for (int x = 0; x != MWIDTH_2; ++x)
                total += r[x];
            for (int x = 0; x != w; ++x) {
                total -= r[x - MWIDTH_2 - 1];
                total += r[x + MWIDTH_2];
                fw[x] = total ? 255 : 0;
            }
        }
    }

    free(row);
    free(fmoving);
}


void noise_dilate(BYTE *VS_RESTRICT moving, const Parms &p) {
    // dilate noise pixels map (with Parms.dilate as range)
    // code is bit modified code of motion_denoise
    // "moving" is actually "noise"  - I do not rename it for simplicity (and lazy)

    if (p.dilate <= 0)
        return;

    const int w = p.pitch;
    const int h = p.height;
    const int w_8 = w / 8;

    //(only dynamic arrays with non-constant size):
    BYTE *fmoving = (BYTE *)malloc(sizeof(BYTE) * h * w);

    // (only dynamic arrays with non-constant size):
    size_t rowpadded_size = ((w + 64) / 8) * 8 + 8;
    signed char *row = (signed char *)_alloca(sizeof(char) * p.y_next * rowpadded_size);

    // we use p.dilate instead of both mheight/2 and mwidth/2
    const int MHEIGHT_2 = p.dilate;
    const int MWIDTH_2 = p.dilate;
    const int Y_NEXT = p.y_next;

    // for noise map expanding, we will not use erode, but only dilate

    // we do not use erode, but we must initialize  fmoving,
      //simply from  noise array data, but of 1 and 0 (not 255 and 0)
      // code added instead of erode
    for (int y = 0; y < h; y++) {
        BYTE *r2 = moving + w * y; // lines
        BYTE *fw = fmoving + w * y;
        for (int x = 0; x < w; x++) {
            fw[x] = r2[x] ? 1 : 0;
        }
    }

    // Dilate
    for (int i = 0; i < Y_NEXT; i++) {
        memzero(&row[i * rowpadded_size], rowpadded_size);
    }

    for (int y = -(MHEIGHT_2)*Y_NEXT; y < h; ++y) {
        signed char *r = row + (y & (Y_NEXT - 1)) + 32;
        int y1 = y - (MHEIGHT_2)*Y_NEXT - Y_NEXT, y2 = y + (MHEIGHT_2)*Y_NEXT;

        motion_summation(fmoving, r, w, h, y1, y2);

        if (y >= 0) {
            BYTE *fw = moving + w * y;
            int total = 0;
            for (int x = 0; x != MWIDTH_2; ++x)
                total += r[x];
            for (int x = 0; x != w; ++x) {
                total -= r[x - MWIDTH_2 - 1];
                total += r[x + MWIDTH_2];
                fw[x] = total ? 255 : 0;
            }
        }
    }

    free(fmoving);
}

static void median_line(const BYTE *p, const BYTE *n, const BYTE *c, const BYTE *c_motion, BYTE *VS_RESTRICT o, int width) {
    for (int x = 0; x < width; ++x) {
        o[x] = c[x];
        if (c_motion[x])
            continue;
        BYTE mn = std::min(p[x], n[x]);
        BYTE mx = std::max(p[x], n[x]);
        o[x] = std::max(mn, o[x]);
        o[x] = std::min(mx, o[x]);
    }
}

void cond_median(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch,
    BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion,
    const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch,
    const Parms &parms) {
        {
            int width = parms.pitch;
            for (int y = 0; y != parms.height; ++y) {
                median_line(p, n, c, c_motion, o, width);

                p += Ppitch;
                c += Cpitch;
                n += Npitch;
                o += Opitch;
                c_noise += parms.pitch;
                c_motion += parms.pitch;
            }
        }
}

static void getminmax(const unsigned char *p, const unsigned char *n,
    unsigned char *VS_RESTRICT minpn, unsigned char *VS_RESTRICT maxpn, int wmod8) {

    if (minpn != 0 && maxpn != 0) {
        for (int x = 0; x < wmod8; ++x) {
            minpn[x] = std::min(p[x], n[x]);  // min from 2 pixels
            maxpn[x] = std::max(p[x], n[x]);  // max from 2 pixels
        }
    } else if (minpn != 0) {
        for (int x = 0; x < wmod8; ++x) {
            minpn[x] = std::min(p[x], n[x]);  // min from 2 pixels
        }
    } else if (maxpn != 0) {
        for (int x = 0; x < wmod8; ++x) {
            maxpn[x] = std::max(p[x], n[x]);  // min from 2 pixels
        }
    }

}

static void find_sdi_p1(const unsigned char *p, const unsigned char *n, const unsigned char *c,
    unsigned char *VS_RESTRICT o, int wmod8) {

    for (int x = 0; x < wmod8; ++x) {
        BYTE minpn = std::min(p[x], n[x]);  // min from 2 pixels
        if (c[x] < minpn)
            o[x] = minpn - c[x];
    }

}


static void find_sdi_m1(const unsigned char *p, const unsigned char *n, const unsigned char *c,
    unsigned char *VS_RESTRICT o, int wmod8) {

    for (int x = 0; x < wmod8; ++x) {
        BYTE maxpn = std::max(p[x], n[x]);  // min from 2 pixels
        if (c[x] > maxpn)
            o[x] = c[x] - maxpn;
    }

}


static void find_sdi_0(const unsigned char *p, const unsigned char *n, const unsigned char *c,
    unsigned char *VS_RESTRICT o, int wmod8) {

    for (int x = 0; x < wmod8; ++x) {
        BYTE minpn = std::min(p[x], n[x]);  // min from 2 pixels
        BYTE maxpn = std::max(p[x], n[x]);  // min from 2 pixels
        if (c[x] < minpn)
            o[x] = minpn - c[x];
        else if (c[x] > maxpn)
            o[x] = c[x] - maxpn;
    }

}

//
// Noise Detection
//

void find_outliers(const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    const BYTE *n, int Npitch, BYTE *VS_RESTRICT o, const Parms &parms) {
    // output o is min difference of cur from prev and next if they both less or both bigger
     // Sign of spot may also be defined
    int width = parms.width;
    int height = parms.height;
    BYTE mn, mx;
    int wmod8 = parms.pitch; // mearest mod8 not less than width and not bigger than pitch
    BYTE *minpn = new BYTE[wmod8 + 16]; // temporary array min(prev[x], next[x]),
    BYTE *maxpn = new BYTE[wmod8 + 16]; // temporary array max(prev[x], next[x]),

    int sign = parms.sign;
    BYTE p1 = parms.p1;

    memzero(o, parms.size);

    for (int y = 0; y != parms.height; ++y) {

        if (parms.ranked) {


            // Simplified Ranked ordered difference detector (S-ROD) for stability to noise
            // Reference:Restoration and Storage of Film and Video Archive Material.
            // P.M.B. van Roosmalen, J. Biemond, and R.L. Lagendijk. 1998.
            // File from Internet: Nato Summer School 1998.pdf
            //
            // Ranked dif with 6 pixels: 3 in prev , 3 in next
            // My modification of article method:
            // I use pixels in same horizontal line, instead of in same vertical row
            // (better for interlaced or telecined video)



                // Set some black-white mode of spots search, from sign parameter:
            switch (sign) {
            case 1://  sign=1 - only black (dark) spots and outliers

                getminmax(p, n, minpn, 0, wmod8);

                for (int x = 1; x < width - 1; ++x) { // omit leftmost and rightmost rows
                    // get extremums of coincident pixels in neighbors frames
                    mn = std::min(minpn[x - 1], minpn[x]);
                    mn = std::min(mn, minpn[x + 1]); // have min from 2+4=6 pixels

                    if (mn > c[x])	o[x] = mn - c[x];
                }
                break;
            case -1://  sign=-1 - only white (light) spots and outliers

                getminmax(p, n, 0, maxpn, wmod8);

                for (int x = 1; x < width - 1; ++x) { // omit leftmost and rightmost rows
                    // get extremums of coincident pixels in neighbors frames

                    mx = std::max(maxpn[x - 1], maxpn[x]);
                    mx = std::max(mx, maxpn[x + 1]);// have max from 2+4=6 pixels


                    if (c[x] > mx)
                        o[x] = c[x] - mx;
                }
                break;
            case 2://  sign=2 - only black (dark) spots, any outliers

                getminmax(p, n, minpn, maxpn, wmod8);

                for (int x = 1; x < width - 1; ++x) { // omit leftmost and rightmost rows
                    // get extremums of coincident pixels in neighbors frames

                    mn = std::min(minpn[x - 1], minpn[x]);
                    mn = std::min(mn, minpn[x + 1]); // have min from 2+4=6 pixels

                    mx = std::max(maxpn[x - 1], maxpn[x]);
                    mx = std::max(mx, maxpn[x + 1]);// have max from 2+4=6 pixels

                    if (c[x] < mn)
                        o[x] = mn - c[x];
                    else if (c[x] > mx)
                        o[x] = std::min<int>(c[x] - mx, p1);
                }
                break;
            case -2:	//  sign=-2 - only white (light) spots, any outliers

                getminmax(p, n, minpn, maxpn, wmod8);

                for (int x = 1; x < width - 1; ++x) { // omit leftmost and rightmost rows
                    // get extremums of coincident pixels in neighbors frames
//
                    mn = std::min(minpn[x - 1], minpn[x]);
                    mn = std::min(mn, minpn[x + 1]); // have min from 2+4=6 pixels

                    mx = std::max(maxpn[x - 1], maxpn[x]);
                    mx = std::max(mx, maxpn[x + 1]);// have max from 2+4=6 pixels


                    if (c[x] > mx)
                        o[x] = c[x] - mx;
                    else if (c[x] < mn)
                        o[x] = std::min<int>(mn - c[x], p1); //3.2
                }
                break;
            case 0:
            default:  //  sign=0 - any spots and outliers (default)
                //  (as used in v. 0.93 code by Kevin Atkinson for all cases)
            // replaced to more fast method (with minmax array) in v.3.2

                getminmax(p, n, minpn, maxpn, wmod8);

                for (int x = 1; x < width - 1; ++x)  // omit leftmost and rightmost rows
                {// get extremums of coincident pixels in neighbors frames

                    mn = std::min(minpn[x - 1], minpn[x]);
                    mn = std::min(mn, minpn[x + 1]); // have min from 2+4=6 pixels

                    mx = std::max(maxpn[x - 1], maxpn[x]);
                    mx = std::max(mx, maxpn[x + 1]);// have max from 2+4=6 pixels

                    if (mn > c[x])	o[x] = mn - c[x];
                    else {
                        if (c[x] > mx)	o[x] = c[x] - mx;
                    }
                }
            }
        }

        else { // old SDI method by 2 pixels
            // Set some black-white mode of spots search, from sign parameter:
            switch (parms.sign) {
            case 1://  sign=1 - only black (dark) spots and outliers

                find_sdi_p1(p, n, c, o, wmod8);

                break;
            case -1://  sign=-1 - only white (light) spots and outliers

                find_sdi_m1(p, n, c, o, wmod8);
                break;
            case 2://  sign=2 - only black (dark) spots, any outliers

                getminmax(p, n, minpn, maxpn, wmod8);

                for (int x = 0; x < width; ++x) {
                    // get extremums of coincident pixels in neighbors frames
                    if (c[x] < minpn[x])
                        o[x] = minpn[x] - c[x];
                    else if (c[x] > maxpn[x])
                        o[x] = std::min<int>(c[x] - maxpn[x], p1);
                }
                break;
            case -2:	//  sign=-2 - only white (light) spots, any outliers

                getminmax(p, n, minpn, maxpn, wmod8);

                for (int x = 0; x < width; ++x) {
                    // get extremums of coincident pixels in neighbors frames
                    if (c[x] > maxpn[x])
                        o[x] = c[x] - maxpn[x];
                    else if (c[x] < minpn[x])
                        o[x] = std::min<int>(minpn[x] - c[x], p1);
                }
                break;
            case 0:
            default:  //  sign=0 - any spots and outliers (default)

                find_sdi_0(p, n, c, o, wmod8);
            }
        }

        p += Ppitch;
        c += Cpitch;
        n += Npitch;
        o += parms.pitch;
    }

    delete[] minpn;
    delete[] maxpn;
}

// change motion mark to 1 if noise
void noise_to_one(BYTE *VS_RESTRICT p_noise, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion, const Parms &parms) {
    int width = parms.width;
    for (int y = 0; y != parms.height; ++y) {
        for (int x = 0; x != width; ++x) {
            if (c_motion[x] && (p_noise[x] || c_noise[x])) c_motion[x] = 1;
            //  motion and  prev or current noise = 1, else no change
        }
        p_noise += parms.pitch;
        c_noise += parms.pitch;
        c_motion += parms.pitch;
    }
}

static void lineblur(const BYTE *c_noise, const BYTE *c_motion, BYTE *VS_RESTRICT o, int width, const Parms &parms) {
    BYTE b[10];	           // buffer
    int xleft = 2 * parms.blur;// sub blur diameter to not exceed bounds
    int xright = width - 2 * parms.blur;

    switch (parms.blur) {
    case 0:                       // no blur
        break;
    case 1:
        for (int x = xleft; x < xright; ++x) {   // sub blur diameter to not exceed bounds
            if (c_noise[x] && !c_motion[x]) {
                if (!c_noise[x - 1]) {            	// spot begins at x
                    b[0] = (o[x - 2] + 2 * o[x - 1] + o[x]) / 4;  // make line (1D) blur with centers at x-1
                    b[1] = (o[x - 1] + 2 * o[x] + o[x + 1]) / 4;	//and at x
                    o[x - 1] = b[0];
                    o[x] = b[1];
                } else if (!c_noise[x + 1]) {          // spot ends at x
                    b[0] = (o[x - 1] + 2 * o[x] + o[x + 1]) / 4;
                    b[1] = (o[x] + 2 * o[x + 1] + o[x + 2]) / 4;
                    o[x] = b[0];
                    o[x + 1] = b[1];
                }
            }
        }
        break;
    case 2:
        for (int x = xleft; x < xright; ++x) {   // sub blur diameter to not exceed bounds
            if (c_noise[x] && !c_motion[x]) {
                if (!c_noise[x - 1]) {            	// spot begins at x
                    b[3] = (o[x - 4] + 2 * o[x - 3] + 3 * o[x - 2] + 2 * o[x - 1] + o[x]) / 9;
                    b[4] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 9;
                    b[5] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 9;
                    b[6] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 9;
                    o[x - 2] = b[3];
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                } else if (!c_noise[x + 1]) {          // spot ends at x
                    b[4] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 9;
                    b[5] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 9;
                    b[6] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 9;
                    b[7] = (o[x] + 2 * o[x + 1] + 3 * o[x + 2] + 2 * o[x + 3] + o[x + 4]) / 9;
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                    o[x + 2] = b[7];
                }
            }
        }
        break;
    case 3:
        for (int x = xleft; x < xright; ++x) {   // sub blur diameter to not exceed bounds
            if (c_noise[x] && !c_motion[x]) {
                if (!c_noise[x - 1]) {            	// spot begins at x
                    b[2] = (o[x - 6] + 2 * o[x - 5] + 3 * o[x - 4] + 4 * o[x - 3] + 3 * o[x - 2] + 2 * o[x - 1] + o[x]) / 16;
                    b[3] = (o[x - 5] + 2 * o[x - 4] + 3 * o[x - 3] + 4 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 16;
                    b[4] = (o[x - 4] + 2 * o[x - 3] + 3 * o[x - 2] + 4 * o[x - 1] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 16;
                    b[5] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 4 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 16;
                    b[6] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 4 * o[x + 1] + 3 * o[x + 2] + 2 * o[x + 3] + o[x + 4]) / 16;
                    b[7] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 4 * o[x + 2] + 3 * o[x + 3] + 2 * o[x + 4] + o[x + 5]) / 16;
                    o[x - 3] = b[2];
                    o[x - 2] = b[3];
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                    o[x + 2] = b[7];
                } else if (!c_noise[x + 1]) {          // spot ends at x
                    b[3] = (o[x - 5] + 2 * o[x - 4] + 3 * o[x - 3] + 4 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 16;
                    b[4] = (o[x - 4] + 2 * o[x - 3] + 3 * o[x - 2] + 4 * o[x - 1] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 16;
                    b[5] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 4 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 16;
                    b[6] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 4 * o[x + 1] + 3 * o[x + 2] + 2 * o[x + 3] + o[x + 4]) / 16;
                    b[7] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 4 * o[x + 2] + 3 * o[x + 3] + 2 * o[x + 4] + o[x + 5]) / 16;
                    b[8] = (o[x] + 2 * o[x + 1] + 3 * o[x + 2] + 4 * o[x + 3] + 3 * o[x + 4] + 2 * o[x + 5] + o[x + 6]) / 16;
                    o[x - 2] = b[3];
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                    o[x + 2] = b[7];
                    o[x + 3] = b[8];
                }
            }
        }
        break;
    case 4:
    default: // limit blur to 4
        for (int x = 8; x < width - 8; ++x) {   // sub blur diameter to not exceed bounds
            if (c_noise[x] && !c_motion[x]) {
                if (!c_noise[x - 1]) {            	// spot begins at x
                    b[1] = (o[x - 8] + 2 * o[x - 7] + 3 * o[x - 6] + 4 * o[x - 5] + 5 * o[x - 4] + 4 * o[x - 3] + 3 * o[x - 2] + 2 * o[x - 1] + o[x]) / 25;
                    b[2] = (o[x - 7] + 2 * o[x - 6] + 3 * o[x - 5] + 4 * o[x - 4] + 5 * o[x - 3] + 4 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 25;
                    b[3] = (o[x - 6] + 2 * o[x - 5] + 3 * o[x - 4] + 4 * o[x - 3] + 5 * o[x - 2] + 4 * o[x - 11] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 25;
                    b[4] = (o[x - 5] + 2 * o[x - 4] + 3 * o[x - 3] + 4 * o[x - 2] + 5 * o[x - 1] + 4 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 25;
                    b[5] = (o[x - 4] + 2 * o[x - 3] + 3 * o[x - 2] + 4 * o[x - 1] + 5 * o[x] + 4 * o[x + 1] + 3 * o[x + 2] + 2 * o[x + 3] + o[x + 4]) / 25;
                    b[6] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 4 * o[x] + 5 * o[x + 1] + 4 * o[x + 2] + 3 * o[x + 3] + 2 * o[x + 4] + o[x + 5]) / 25;
                    b[7] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 4 * o[x + 1] + 5 * o[x + 2] + 4 * o[x + 3] + 3 * o[x + 4] + 2 * o[x + 5] + o[x + 6]) / 25;
                    b[8] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 4 * o[x + 2] + 5 * o[x + 3] + 4 * o[x + 4] + 3 * o[x + 5] + 2 * o[x + 6] + o[x + 7]) / 25;
                    b[9] = (o[x] + 2 * o[x + 1] + 3 * o[x + 2] + 4 * o[x + 3] + 5 * o[x + 4] + 4 * o[x + 5] + 3 * o[x + 6] + 2 * o[x + 7] + o[x + 8]) / 25;
                    o[x - 4] = b[1];
                    o[x - 3] = b[2];
                    o[x - 2] = b[3];
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                    o[x + 2] = b[7];
                    o[x + 3] = b[8];
                } else if (!c_noise[x + 1]) {          // spot ends at x
                    b[2] = (o[x - 7] + 2 * o[x - 6] + 3 * o[x - 5] + 4 * o[x - 4] + 5 * o[x - 3] + 4 * o[x - 2] + 3 * o[x - 1] + 2 * o[x] + o[x + 1]) / 25;
                    b[3] = (o[x - 6] + 2 * o[x - 5] + 3 * o[x - 4] + 4 * o[x - 3] + 5 * o[x - 2] + 4 * o[x - 11] + 3 * o[x] + 2 * o[x + 1] + o[x + 2]) / 25;
                    b[4] = (o[x - 5] + 2 * o[x - 4] + 3 * o[x - 3] + 4 * o[x - 2] + 5 * o[x - 1] + 4 * o[x] + 3 * o[x + 1] + 2 * o[x + 2] + o[x + 3]) / 25;
                    b[5] = (o[x - 4] + 2 * o[x - 3] + 3 * o[x - 2] + 4 * o[x - 1] + 5 * o[x] + 4 * o[x + 1] + 3 * o[x + 2] + 2 * o[x + 3] + o[x + 4]) / 25;
                    b[6] = (o[x - 3] + 2 * o[x - 2] + 3 * o[x - 1] + 4 * o[x] + 5 * o[x + 1] + 4 * o[x + 2] + 3 * o[x + 3] + 2 * o[x + 4] + o[x + 5]) / 25;
                    b[7] = (o[x - 2] + 2 * o[x - 1] + 3 * o[x] + 4 * o[x + 1] + 5 * o[x + 2] + 4 * o[x + 3] + 3 * o[x + 4] + 2 * o[x + 5] + o[x + 6]) / 25;
                    b[8] = (o[x - 1] + 2 * o[x] + 3 * o[x + 1] + 4 * o[x + 2] + 5 * o[x + 3] + 4 * o[x + 4] + 3 * o[x + 5] + 2 * o[x + 6] + o[x + 7]) / 25;
                    b[9] = (o[x] + 2 * o[x + 1] + 3 * o[x + 2] + 4 * o[x + 3] + 5 * o[x + 4] + 4 * o[x + 5] + 3 * o[x + 6] + 2 * o[x + 7] + o[x + 8]) / 25;
                    o[x - 3] = b[2];
                    o[x - 2] = b[3];
                    o[x - 1] = b[4];
                    o[x] = b[5];
                    o[x + 1] = b[6];
                    o[x + 2] = b[7];
                    o[x + 3] = b[8];
                    o[x + 4] = b[9];
                }
            }
        }
    }

}

static void tsmooth_line(const BYTE *p, const BYTE *n, BYTE *VS_RESTRICT o, const BYTE *c_noise, const BYTE *c_motion, int width, int tsmoothLimit) {

    for (int x = 0; x != width; ++x) {
        if (!c_noise[x] && !c_motion[x]) {
            // temporal smoothing of almost static pixels (no noise(spots), no motion).
            // Algo is similar to FaeryDust?
            int m = (p[x] + o[x] + n[x]) / 3;        // mean value (integer precise)
            int corr = m - o[x];
            if (corr > (tsmoothLimit * 2)) corr = 0;
            if (corr < -(tsmoothLimit * 2)) corr = 0;
            if (corr > tsmoothLimit) corr = (tsmoothLimit * 2) - corr;
            if (corr < -tsmoothLimit) corr = -(tsmoothLimit * 2) - corr;
            o[x] += corr;
        }
    }
}

void remove_outliers(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch,
    BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion,
    const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch,
    const Parms &parms) {
        {
            int width = parms.width;
            int tsmooth = parms.tsmooth;

            for (int y = 0; y != parms.height; ++y) {

                if (parms.fitluma) {
                    int lumac = 0;
                    int lumad = 0;
                    for (int x = 0; x != width; ++x) {
                        lumac += c[x];               //  prepare line lumas
                        lumad += p[x] + n[x];
                    }
                    lumac /= width;
                    lumad /= (width * 2);
                    lumad = lumac - lumad;      //  luma delta

                    for (int x = 0; x != width; ++x) {
                        if (c_noise[x] && !c_motion[x]) {
                            // luma correction  of deleted spot place
                            int ox = ((p[x] + n[x]) >> 1) + lumad;
                            o[x] = std::min(std::max(ox, 0), 255);
                        } else
                            o[x] = c[x];
                    }
                } else // no fitluma
                {
                    for (int x = 0; x != width; ++x) {
                        if (c_noise[x] && !c_motion[x]) {
                            o[x] = (p[x] + n[x]) >> 1;
                        } else
                            o[x] = c[x];
                    }

                }

                // blur line near spot segments
                if (parms.blur) lineblur(c_noise, c_motion, o, width, parms);

                if (parms.tsmooth) { // temporal smoothing
                    // must be after segments blur!
                    tsmooth_line(p, n, o, c_noise, c_motion, width, tsmooth);
                }

                p += Ppitch;
                c += Cpitch;
                n += Npitch;
                o += Opitch;
                c_noise += parms.pitch;
                c_motion += parms.pitch;
            }

        }
}

void mark_outliers(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion,
    const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch,
    const Parms &parms) {
    for (int y = 0; y != parms.height; ++y) {
        for (int x = 0; x != parms.width; ++x) {
            if (c_noise[x] && !c_motion[x])
                o[x] = parms.mark_v;
            else if (!c_noise[x] && c_motion[x])
                o[x] = 126 + c[x] / 2;
            else if (c[x] == parms.mark_v) {
                if (c[x] >= 1)  o[x] = c[x] - 1;
                else o[x] = c[x] + 1;

            } else
                o[x] = c[x];
        }
        p += Ppitch;
        c += Cpitch;
        n += Npitch;
        o += Opitch;
        c_noise += parms.pitch;
        c_motion += parms.pitch;
    }
}

void map_outliers(const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch,
    BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion,
    const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch,
    const Parms &parms) {
    for (int y = 0; y != parms.height; ++y) {
        for (int x = 0; x != parms.width; ++x) {
            if ((c_noise != 0) && c_noise[x]) {
                if (c_motion[x])
                    o[x] = 159;
                else
                    o[x] = 255;
            } else if (c_motion[x])
                o[x] = 95;
            else
                o[x] = 0;

        }
        p += Ppitch;
        c += Cpitch;
        n += Npitch;
        o += Opitch;
        if (c_noise)
            c_noise += parms.pitch;
        c_motion += parms.pitch;
    }
}

void mark_segments(Segments &segs, const BYTE *p, int Ppitch,
    const BYTE *c, int Cpitch, BYTE *c_noise, BYTE *c_motion,
    const BYTE *n, int Npitch,
    BYTE *o, int Opitch,
    const Parms &parms) {

    for (int y = 0; y != parms.height; ++y) {
        for (int x = 0; x != parms.width; ++x) {
            if (c_noise[x] && !c_motion[x])
                o[x] = parms.mark_v;
            else if (!c_noise[x] && c_motion[x])
                o[x] = 126 + c[x] / 2;
            else if (c[x] == parms.mark_v) {
                if (c[x] >= 1)  o[x] = c[x] - 1;
                else o[x] = c[x] + 1;

            } else  o[x] = c[x];
        }
        p += Ppitch;
        c += Cpitch;
        n += Npitch;
        o += Opitch;
        c_noise += parms.pitch;
        c_motion += parms.pitch;
    }
}

void motion_merge(BYTE *VS_RESTRICT c_motion, BYTE *VS_RESTRICT n_motion, BYTE *VS_RESTRICT m_motion, const Parms &parms) {
    for (int y = 0; y != parms.height; ++y) {
        for (int x = 0; x != parms.width; ++x) {
            m_motion[x] = (c_motion[x] || n_motion[x]) ? 255 : 0;
        }

        c_motion += parms.pitch;
        n_motion += parms.pitch;
        m_motion += parms.pitch;
    }
}


//    1D  line triangle blur near deleted spot segments.

static void segment_blur(BYTE *VS_RESTRICT o, int x1, int x2, int blur, BYTE *VS_RESTRICT b) {
    //    int blur = parms.blur; // 1d blur radius
    //    b must be BYTE buffer with size >=16

    switch (blur) {
    case 0:                       // no blur
        break;
    case 1:
        b[0] = (o[x1 - 2] + 2 * o[x1 - 1] + o[x1]) / 4;  // make line (1D) blur with centers at x-1
        b[1] = (o[x1 - 1] + 2 * o[x1] + o[x1 + 1]) / 4;	//and at x
        o[x1 - 1] = b[0];
        o[x1] = b[1];
        b[2] = (o[x2 - 1] + 2 * o[x2] + o[x2 + 1]) / 4;
        b[3] = (o[x2] + 2 * o[x2 + 1] + o[x2 + 2]) / 4;
        o[x2] = b[2];
        o[x2 + 1] = b[3];
        break;
    case 2:
        b[0] = (o[x1 - 4] + 2 * o[x1 - 3] + 3 * o[x1 - 2] + 2 * o[x1 - 1] + o[x1]) / 9;
        b[1] = (o[x1 - 3] + 2 * o[x1 - 2] + 3 * o[x1 - 1] + 2 * o[x1] + o[x1 + 1]) / 9;
        b[2] = (o[x1 - 2] + 2 * o[x1 - 1] + 3 * o[x1] + 2 * o[x1 + 1] + o[x1 + 2]) / 9;
        b[3] = (o[x1 - 1] + 2 * o[x1] + 3 * o[x1 + 1] + 2 * o[x1 + 2] + o[x1 + 3]) / 9;
        o[x1 - 2] = b[0];
        o[x1 - 1] = b[1];
        o[x1] = b[2];
        o[x1 + 1] = b[3];
        b[4] = (o[x2 - 3] + 2 * o[x2 - 2] + 3 * o[x2 - 1] + 2 * o[x2] + o[x2 + 1]) / 9;
        b[5] = (o[x2 - 2] + 2 * o[x2 - 1] + 3 * o[x2] + 2 * o[x2 + 1] + o[x2 + 2]) / 9;
        b[6] = (o[x2 - 1] + 2 * o[x2] + 3 * o[x2 + 1] + 2 * o[x2 + 2] + o[x2 + 3]) / 9;
        b[7] = (o[x2] + 2 * o[x2 + 1] + 3 * o[x2 + 2] + 2 * o[x2 + 3] + o[x2 + 4]) / 9;
        o[x2 - 1] = b[4];
        o[x2] = b[5];
        o[x2 + 1] = b[6];
        o[x2 + 2] = b[7];
        break;
    case 3:
        b[0] = (o[x1 - 6] + 2 * o[x1 - 5] + 3 * o[x1 - 4] + 4 * o[x1 - 3] + 3 * o[x1 - 2] + 2 * o[x1 - 1] + o[x1]) / 16;
        b[1] = (o[x1 - 5] + 2 * o[x1 - 4] + 3 * o[x1 - 3] + 4 * o[x1 - 2] + 3 * o[x1 - 1] + 2 * o[x1] + o[x1 + 1]) / 16;
        b[2] = (o[x1 - 4] + 2 * o[x1 - 3] + 3 * o[x1 - 2] + 4 * o[x1 - 1] + 3 * o[x1] + 2 * o[x1 + 1] + o[x1 + 2]) / 16;
        b[3] = (o[x1 - 3] + 2 * o[x1 - 2] + 3 * o[x1 - 1] + 4 * o[x1] + 3 * o[x1 + 1] + 2 * o[x1 + 2] + o[x1 + 3]) / 16;
        b[4] = (o[x1 - 2] + 2 * o[x1 - 1] + 3 * o[x1] + 4 * o[x1 + 1] + 3 * o[x1 + 2] + 2 * o[x1 + 3] + o[x1 + 4]) / 16;
        b[5] = (o[x1 - 1] + 2 * o[x1] + 3 * o[x1 + 1] + 4 * o[x1 + 2] + 3 * o[x1 + 3] + 2 * o[x1 + 4] + o[x1 + 5]) / 16;
        o[x1 - 3] = b[0];
        o[x1 - 2] = b[1];
        o[x1 - 1] = b[2];
        o[x1] = b[3];
        o[x1 + 1] = b[4];
        o[x1 + 2] = b[5];
        b[6] = (o[x2 - 5] + 2 * o[x2 - 4] + 3 * o[x2 - 3] + 4 * o[x2 - 2] + 3 * o[x2 - 1] + 2 * o[x2] + o[x2 + 1]) / 16;
        b[7] = (o[x2 - 4] + 2 * o[x2 - 3] + 3 * o[x2 - 2] + 4 * o[x2 - 1] + 3 * o[x2] + 2 * o[x2 + 1] + o[x2 + 2]) / 16;
        b[8] = (o[x2 - 3] + 2 * o[x2 - 2] + 3 * o[x2 - 1] + 4 * o[x2] + 3 * o[x2 + 1] + 2 * o[x2 + 2] + o[x2 + 3]) / 16;
        b[9] = (o[x2 - 2] + 2 * o[x2 - 1] + 3 * o[x2] + 4 * o[x2 + 1] + 3 * o[x2 + 2] + 2 * o[x2 + 3] + o[x2 + 4]) / 16;
        b[10] = (o[x2 - 1] + 2 * o[x2] + 3 * o[x2 + 1] + 4 * o[x2 + 2] + 3 * o[x2 + 3] + 2 * o[x2 + 4] + o[x2 + 5]) / 16;
        b[11] = (o[x2] + 2 * o[x2 + 1] + 3 * o[x2 + 2] + 4 * o[x2 + 3] + 3 * o[x2 + 4] + 2 * o[x2 + 5] + o[x2 + 6]) / 16;
        o[x2 - 2] = b[6];
        o[x2 - 1] = b[7];
        o[x2] = b[8];
        o[x2 + 1] = b[9];
        o[x2 + 2] = b[10];
        o[x2 + 3] = b[11];
        break;
    case 4:
    default: // limit blur range to 4
        b[0] = (o[x1 - 8] + 2 * o[x1 - 7] + 3 * o[x1 - 6] + 4 * o[x1 - 5] + 5 * o[x1 - 4] + 4 * o[x1 - 3] + 3 * o[x1 - 2] + 2 * o[x1 - 1] + o[x1]) / 25;
        b[1] = (o[x1 - 7] + 2 * o[x1 - 6] + 3 * o[x1 - 5] + 4 * o[x1 - 4] + 5 * o[x1 - 3] + 4 * o[x1 - 2] + 3 * o[x1 - 1] + 2 * o[x1] + o[x1 + 1]) / 25;
        b[2] = (o[x1 - 6] + 2 * o[x1 - 5] + 3 * o[x1 - 4] + 4 * o[x1 - 3] + 5 * o[x1 - 2] + 4 * o[x1 - 11] + 3 * o[x1] + 2 * o[x1 + 1] + o[x1 + 2]) / 25;
        b[3] = (o[x1 - 5] + 2 * o[x1 - 4] + 3 * o[x1 - 3] + 4 * o[x1 - 2] + 5 * o[x1 - 1] + 4 * o[x1] + 3 * o[x1 + 1] + 2 * o[x1 + 2] + o[x1 + 3]) / 25;
        b[4] = (o[x1 - 4] + 2 * o[x1 - 3] + 3 * o[x1 - 2] + 4 * o[x1 - 1] + 5 * o[x1] + 4 * o[x1 + 1] + 3 * o[x1 + 2] + 2 * o[x1 + 3] + o[x1 + 4]) / 25;
        b[5] = (o[x1 - 3] + 2 * o[x1 - 2] + 3 * o[x1 - 1] + 4 * o[x1] + 5 * o[x1 + 1] + 4 * o[x1 + 2] + 3 * o[x1 + 3] + 2 * o[x1 + 4] + o[x1 + 5]) / 25;
        b[6] = (o[x1 - 2] + 2 * o[x1 - 1] + 3 * o[x1] + 4 * o[x1 + 1] + 5 * o[x1 + 2] + 4 * o[x1 + 3] + 3 * o[x1 + 4] + 2 * o[x1 + 5] + o[x1 + 6]) / 25;
        b[7] = (o[x1 - 1] + 2 * o[x1] + 3 * o[x1 + 1] + 4 * o[x1 + 2] + 5 * o[x1 + 3] + 4 * o[x1 + 4] + 3 * o[x1 + 5] + 2 * o[x1 + 6] + o[x1 + 7]) / 25;
        b[8] = (o[x1] + 2 * o[x1 + 1] + 3 * o[x1 + 2] + 4 * o[x1 + 3] + 5 * o[x1 + 4] + 4 * o[x1 + 5] + 3 * o[x1 + 6] + 2 * o[x1 + 7] + o[x1 + 8]) / 25;
        o[x1 - 4] = b[0];
        o[x1 - 3] = b[1];
        o[x1 - 2] = b[2];
        o[x1 - 1] = b[3];
        o[x1] = b[4];
        o[x1 + 1] = b[5];
        o[x1 + 2] = b[6];
        o[x1 + 3] = b[7];
        b[8] = (o[x2 - 7] + 2 * o[x2 - 6] + 3 * o[x2 - 5] + 4 * o[x2 - 4] + 5 * o[x2 - 3] + 4 * o[x2 - 2] + 3 * o[x2 - 1] + 2 * o[x2] + o[x2 + 1]) / 25;
        b[9] = (o[x2 - 6] + 2 * o[x2 - 5] + 3 * o[x2 - 4] + 4 * o[x2 - 3] + 5 * o[x2 - 2] + 4 * o[x2 - 11] + 3 * o[x2] + 2 * o[x2 + 1] + o[x2 + 2]) / 25;
        b[10] = (o[x2 - 5] + 2 * o[x2 - 4] + 3 * o[x2 - 3] + 4 * o[x2 - 2] + 5 * o[x2 - 1] + 4 * o[x2] + 3 * o[x2 + 1] + 2 * o[x2 + 2] + o[x2 + 3]) / 25;
        b[11] = (o[x2 - 4] + 2 * o[x2 - 3] + 3 * o[x2 - 2] + 4 * o[x2 - 1] + 5 * o[x2] + 4 * o[x2 + 1] + 3 * o[x2 + 2] + 2 * o[x2 + 3] + o[x2 + 4]) / 25;
        b[12] = (o[x2 - 3] + 2 * o[x2 - 2] + 3 * o[x2 - 1] + 4 * o[x2] + 5 * o[x2 + 1] + 4 * o[x2 + 2] + 3 * o[x2 + 3] + 2 * o[x2 + 4] + o[x2 + 5]) / 25;
        b[13] = (o[x2 - 2] + 2 * o[x2 - 1] + 3 * o[x2] + 4 * o[x2 + 1] + 5 * o[x2 + 2] + 4 * o[x2 + 3] + 3 * o[x2 + 4] + 2 * o[x2 + 5] + o[x2 + 6]) / 25;
        b[14] = (o[x2 - 1] + 2 * o[x2] + 3 * o[x2 + 1] + 4 * o[x2 + 2] + 5 * o[x2 + 3] + 4 * o[x2 + 4] + 3 * o[x2 + 5] + 2 * o[x2 + 6] + o[x2 + 7]) / 25;
        b[15] = (o[x2] + 2 * o[x2 + 1] + 3 * o[x2 + 2] + 4 * o[x2 + 3] + 5 * o[x2 + 4] + 4 * o[x2 + 5] + 3 * o[x2 + 6] + 2 * o[x2 + 7] + o[x2 + 8]) / 25;
        o[x2 - 3] = b[8];
        o[x2 - 2] = b[9];
        o[x2 - 1] = b[10];
        o[x2] = b[11];
        o[x2 + 1] = b[12];
        o[x2 + 2] = b[13];
        o[x2 + 3] = b[14];
        o[x2 + 4] = b[15];
    }

}


//
// temporal smoothing

void temporal_smooth(Segments &segs, const BYTE *p, int Ppitch, const BYTE *n, int Npitch,
    BYTE *VS_RESTRICT o, int Opitch, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT c_motion, const Parms &parms) {
        {
            int width = parms.width;
            int sigmann2 = (parms.tsmooth) * (parms.tsmooth); // approximate value of noise variance for Wiener temporal smoothing
            int tsmooth = parms.tsmooth;

            Segment *seg = segs.data;
            segs.end->ly = parms.height; // to avoid having to check for < segs.end

            for (int y = 0; y < parms.height; ++y) {
                for (; seg->ly == y; ++seg) {

                    int x1 = std::max<int>(seg->lx1, 0);
                    int x2 = std::min<int>(seg->lx2, width - 1);
                    for (int x = x1; x <= x2; ++x) {
                        c_noise[x] = 255;
                    }
                }

                tsmooth_line(p, n, o, c_noise, c_motion, width, tsmooth);

                p += Ppitch;
                n += Npitch;
                o += Opitch;
                c_noise += parms.pitch;
                c_motion += parms.pitch;
            }

        }
}

void remove_segments(Segments &segs, const BYTE *p, int Ppitch, const BYTE *c, int Cpitch,
    const BYTE *n, int Npitch, BYTE *VS_RESTRICT o, int Opitch, BYTE *VS_RESTRICT c_noise, const Parms &parms) {

    int width = parms.width;
    int blur = parms.blur;
    BYTE blurBuffer[16]; // buffer for blur

    bool fitluma = parms.fitluma;
    int lumaAdd;
    int dilateX = abs(parms.dilate);
    int dilateY = abs(parms.dilate / parms.y_next);

    int dy;
    Segment *seg = segs.data;
    segs.end->ly = parms.height; // to avoid having to check for < segs.end

    // for adaptive dilate value
    int dilateX1, dilateX2, dilateYC;
    int similar = parms.p2;

    int opitch = Opitch * parms.y_next;
    int ppitch = Ppitch * parms.y_next;
    int npitch = Npitch * parms.y_next;
    int cnpitch = parms.pitch * parms.y_next; // noise

    memzero(c_noise, parms.size);

    o += opitch * dilateY; // skip first dilateY lines
    p += ppitch * dilateY; //
    n += npitch * dilateY; //
    c_noise += cnpitch * dilateY;


    for (int y = dilateY * parms.y_next; y < parms.height - dilateY * parms.y_next; ++y) {
        for (; seg->ly < y; ++seg) {   //skip segnents in skipped dilateY lines
            ;
        }

        for (; seg->ly == y; ++seg) {
            if (seg->data.s.is_noise) {
                int x1 = std::max<int>(seg->lx1, dilateX + 1 + blur);
                int x2 = std::min<int>(seg->lx2, width - dilateX - 2 - blur);
                if (parms.dilate < 0) {
                    // constrained dilation
                    dilateX1 = 0;
                    for (int i = 1; i <= dilateX; i++) {
                        if (abs(o[x1 - i] - o[x1]) > similar) break;
                        dilateX1 = i;
                    }
                    dilateX2 = 0;
                    for (int i = 1; i <= dilateX; i++) {
                        if (abs(o[x1 + i] - o[x1]) > similar) break;
                        dilateX2 = i;
                    }
                    x1 -= dilateX1;
                    x2 += dilateX2;
                    dilateYC = std::min(dilateY, std::max(dilateX1, dilateX2) / parms.y_next);
                } else {
                    // unconstrained dilation
                    x1 -= dilateX;
                    x2 += dilateX;
                    dilateYC = dilateY;
                }

                // local (segment) luma correction to fit luma to neighbor pixels
                if (fitluma) {
                    if (y >= parms.y_next && y < parms.height - parms.y_next) {// middle lines
                        lumaAdd = o[x1 - 1] + o[x2 + 1] - (p[x1 - 1] + p[x2 + 1] + n[x1 - 1] + n[x2 + 1]) +
                            (o[x1 - 1 + opitch] + o[x2 + 1 + opitch] + o[x1 - 1 - opitch] + o[x2 + 1 - opitch]) / 2;
                        lumaAdd /= 4;
                    } else {
                        lumaAdd = o[x1 - 1] + o[x2 + 1] - (p[x1 - 1] + p[x2 + 1] + n[x1 - 1] + n[x2 + 1]) / 2;
                        lumaAdd /= 2;
                    }
                } else
                    lumaAdd = 0;

                for (int x = x1; x <= x2; ++x) {
                    o[x] = std::min(std::max((p[x] + n[x]) / 2 + lumaAdd, 0), 255);
                    c_noise[x] = 255;
                }

                if (parms.blur) segment_blur(o, x1, x2, blur, blurBuffer);

                if (dilateYC > 0) {
                    if (!(seg->what & B_LINK)) {// the most top master (or single), dilate to top
                        for (dy = 0; dy < dilateYC; dy++) {
                            o += (-opitch); // prev line
                            p += (-ppitch);
                            n += (-npitch);
                            c_noise += (-cnpitch);
                            for (int x = x1; x <= x2; ++x) {
                                o[x] = std::min(std::max((p[x] + n[x]) / 2 + lumaAdd, 0), 255);
                                c_noise[x] = 255;
                            }
                            if (parms.blur) segment_blur(o, x1, x2, blur, blurBuffer);
                        }
                        o += opitch * dilateYC; // restore
                        p += ppitch * dilateYC;
                        n += npitch * dilateYC;
                        c_noise += cnpitch * dilateYC;
                    }
                    if (!(seg->what & B_MASTER)) {   // the most bottom link (or single), dilate to bottom
                        for (dy = 0; dy < dilateYC; dy++) {
                            o += opitch; // next line
                            p += ppitch;
                            n += npitch;
                            c_noise += cnpitch;
                            for (int x = x1; x <= x2; ++x) {
                                o[x] = std::min(std::max((p[x] + n[x]) / 2 + lumaAdd, 0), 255);
                                c_noise[x] = 255;
                            }
                            if (parms.blur) segment_blur(o, x1, x2, blur, blurBuffer);
                        }
                        o += -opitch * dilateYC;// restore
                        p += -ppitch * dilateYC;
                        n += -npitch * dilateYC;
                        c_noise += -cnpitch * dilateYC;
                    }
                }
            }
        }
        p += Ppitch;
        c += Cpitch;
        n += Npitch;
        o += Opitch;
        c_noise += parms.pitch;
    }
}


//
// set color at places of deleted spot to mean value
//
void clean_color_plane(const BYTE *p, int ppitch, const BYTE *c, int cpitch, const BYTE *n, int npitch, int widthUV, int heightUV, BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT o, int opitch, Parms &Params)
{
    int mult = (1 << 10) / 3;
    int y_next = Params.y_next;

    if (heightUV == Params.height / 2 && widthUV == Params.width / 2) // YV12
    {
        int noisePitch = Params.pitch * y_next;
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 4 points of luma noise
                if (c_noise[i * 2] || c_noise[i * 2 + 1] || c_noise[i * 2 + noisePitch] || c_noise[i * 2 + noisePitch + 1]) {
                    o[i] = (p[i] + n[i] + c[i]) / 3;
                } else {
                    o[i] = c[i];
                }
            }
            p += ppitch;
            c += cpitch;
            n += npitch;
            o += opitch;
            c_noise += Params.pitch * (1 + (y_next % 2) + 2 * (j % y_next));

        }
    } else if (heightUV == Params.height && widthUV == Params.width / 2) // YV16
    {
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 2 points of luma noise
                if (c_noise[i * 2] || c_noise[i * 2 + 1]) {
                    o[i] = (p[i] + n[i] + c[i]) / 3;
                } else {
                    o[i] = c[i];
                }
            }
            p += ppitch;
            c += cpitch;
            n += npitch;
            o += opitch;
            c_noise += Params.pitch;

        }
    } else if (heightUV == Params.height && widthUV == Params.width) // YV24
    {
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 2 points of luma noise
                if (c_noise[i * 2]) {
                    o[i] = (p[i] + n[i] + c[i]) / 3;
                } else {
                    o[i] = c[i];
                }
            }
            p += ppitch;
            c += cpitch;
            n += npitch;
            o += opitch;
            c_noise += Params.pitch;

        }
    }
}

//
// mark spots by color
//

void mark_color_plane(BYTE *VS_RESTRICT c_noise, BYTE *VS_RESTRICT ptrV, int pitchV, int widthUV, int heightUV, Parms &Params) {
    int mark_v = Params.mark_v;
    int y_next = Params.y_next;
    BYTE mark_color = ((mark_v + 2) % 3) * mark_v / 2;

    if (heightUV == Params.height / 2 && widthUV == Params.width / 2) // YV12
    {
        int noisePitch = Params.pitch * y_next;
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 4 points of luma noise
                if ((c_noise[i * 2] == 255)
                    || (c_noise[i * 2 + 1] == 255)
                    || (c_noise[noisePitch + i * 2] == 255)
                    || (c_noise[noisePitch + i * 2 + 1] == 255)) {
                    ptrV[i] = mark_color; // mark spot by color
                }
            }
            ptrV += pitchV;
            c_noise += Params.pitch * (1 + (y_next % 2) + 2 * (j % y_next));
        }
    } else if (heightUV == Params.height && widthUV == Params.width / 2) // YV16
    {
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 2 points of luma noise
                if ((c_noise[i * 2] == 255)
                    || (c_noise[i * 2 + 1] == 255)) {
                    ptrV[i] = mark_color; // mark spot by color
                }
            }
            ptrV += pitchV;
            c_noise += Params.pitch;
        }
    } else if (heightUV == Params.height && widthUV == Params.width) // YV24
    {
        for (int j = 0; j < heightUV; j++) {
            for (int i = 0; i < widthUV; i++) {
                // check 2 points of luma noise
                if (c_noise[i * 2] == 255) {
                    ptrV[i] = mark_color; // mark spot by color
                }
            }
            ptrV += pitchV;
            c_noise += Params.pitch;
        }
    }
}

static int motion_count(const BYTE *motion, int height, int width, int pitch) {
    // counter of motion points
    int count = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (motion[x]) count++;
        }
        motion += pitch;
    }

    return count;
}

//
// set full motion map to motion at scenechange (if motion points percent > threshold)
//

void motion_scene(BYTE *VS_RESTRICT motion, const Parms &parms) {
    int width_mod8 = (parms.width / 8) * 8;

    int count = motion_count(motion, parms.height, width_mod8, parms.pitch);

    int mpercent = (count * 100) / (parms.height * width_mod8);

    if (mpercent > parms.mscene) {// new scene found, set all points as motion
        size_t fullsize = parms.height * parms.pitch;
        for (int h = 0; h < parms.height; h++) {
            memset(motion, 255, parms.width);
            motion += parms.pitch;
        }
    }
}

void add_external_mask(const BYTE *ext, int ext_pitch, BYTE *VS_RESTRICT motion, int pitch, int width, int height) {	// add external mask of good objects (for example SAD<threshold) to motion
    // motion on input is binary mask 0 or 255
    //  result = ext>127 || motion

    for (int h = 0; h < height; h++) {
        for (int w = 0; w < width; w++) {
            if (ext[w] > 127)
                motion[w] = 255; // set to 255 if motion OR extMask, set to 0 if not
        }
        motion += pitch;
        ext += ext_pitch;

    }
}

static void	print_timestamp(char buf_0[], double t) {
    const int		all = int(floor(t * 100 + 0.5));
    const int		cs = all % 100;
    const int		s = (all / 100) % 60;
    const int		m = (all / 6000) % 60;
    const int		h = all / 360000;

    sprintf(buf_0, "%01d:%02d:%02d.%02d", h, m, s, cs);
}


void print_segments(const Segments &segs, FILE *f_ptr, int fn, double frate, int w, int h, int spotmax1, int spotmax2) {
    // First, checks if this is a cluster of false positive
    int				count = 0;
    double			sum1[2] = { 0, 0 };
    double			sum2[2] = { 0, 0 };
    for (const Segment *cur = segs.data; cur != segs.end; ++cur) {
        if ((cur->what & B_LINK) == 0) {
            const Segment::Data &data = cur->data;
            const double	x = (data.b.x1 + data.b.x2) * 0.5;
            const double	y = (data.b.y1 + data.b.y2) * 0.5;
            sum1[0] += x;
            sum1[1] += y;
            sum2[0] += x * x;
            sum2[1] += y * y;
            ++count;
        }
    }
    const double	std_dev[2] =
    {
        sqrt(sum2[0] - sum1[0] * sum1[0]),
        sqrt(sum2[1] - sum1[1] * sum1[1])
    };
    const double	std_dev_x = std_dev[0] / w;
    const double	std_dev_y = std_dev[1] / h;
    const double	std_dev_diag =
        sqrt((std_dev_x * std_dev_x + std_dev_y * std_dev_y) * 0.5);

    const double	thr = 0.5 / sqrt(12.0);
    if ((count > spotmax1 && std_dev_diag < thr)
        || count > spotmax2) {
        return;
    }
    // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

    const double	tb = std::max((fn - 0.5) / frate, 0.0);
    const double	te = (fn + 0.5) / frate;
    char				tb_0[63 + 1];
    char				te_0[63 + 1];
    print_timestamp(tb_0, tb);
    print_timestamp(te_0, te);

    for (const Segment *cur = segs.data; cur != segs.end; ++cur) {
        if ((cur->what & B_LINK) == 0) {
            const Segment::Data &data = cur->data;
            if (data.s.is_noise) {
                const int		margin = 2;
                const int		left = std::min(std::max(data.b.x1 - margin, 0), w);
                const int		right = std::min(std::max(data.b.x2 + margin, 0), w);
                const int		top = std::min(std::max(data.b.y1 - margin, 0), h);
                const int		bottom = std::min(std::max(data.b.y2 + margin, 0), h);
                const int		cx = (data.b.x1 + data.b.x2 + 1) >> 1;
                const int		cy = (data.b.y1 + data.b.y2 + 1) >> 1;
                fprintf(f_ptr,
                    "Dialogue: 0,%s,%s,Mask,,0,0,0,,"
                    "{\\pos(%d,%d)\\fscx100\\fscy100\\p1}m %d %d l %d %d %d %d %d %d{\\p0}\n",
                    tb_0, te_0,
                    cx, cy,
                    left - cx, top - cy,
                    right - cx, top - cy,
                    right - cx, bottom - cy,
                    left - cx, bottom - cy
                );
            }
        }
    }
}

