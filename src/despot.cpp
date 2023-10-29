/*
    DeSpot - Conditional Temporal Despotting Filter for Avisynth 2.5

    Main file

    Version 3.5.3 October 23, 2008.



Copyright (C)2003-2008 Alexander G. Balakhnin aka Fizick.
bag@hotmail.ru         http://bag.hotmail.ru
under the GNU General Public Licence version 2.

This plugin is based on  Conditional Temporal Median Filter (c-plugin) for Avisynth 2.5
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

#include <avisynth.h>
#include <cstdlib>
#include <cstring>
#include <cmath>

struct WorkingData {
    int num = -1;
    BYTE *noise = nullptr;
    BYTE *motion = nullptr;
    Segment *segments = nullptr;
    void init(size_t);
    ~WorkingData();
};

struct Frame {
    int num = -1;
    WorkingData *data = nullptr;
    PVideoFrame frame;
    const BYTE *y = nullptr; // Y plane of frame
    int pitch = 0; // pitch of Y plane
    BYTE *noise = nullptr;
    BYTE *motion = 0;
    Segments segments;
};

constexpr size_t fcache_size = 16;

struct Buffer {
    PClip childclip;
    WorkingData data_cache[2];
    Frame frame_cache[fcache_size];
    BYTE *motion; // summary motion
    int num_frames;

    Frame *get_frame(int n, IScriptEnvironment *env) {
        if (n < 0 || n >= num_frames)
            return 0;
        Frame *f = frame_cache + n % fcache_size;
        if (f->num != n)
            ready_frame(n, f, env);
        return f;
    }

    void alloc_noise(Frame *f) {
        if (!f->data) ready_data(f);
        f->noise = f->data->noise;
    }

    void alloc_motion(Frame *f) {
        if (!f->data) ready_data(f);
        f->motion = f->data->motion;
    }

    void alloc_segments(Frame *f) {
        if (!f->data) ready_data(f);
        f->segments.data = f->data->segments;
        f->segments.clear();
    }

    void ready_frame(int n, Frame *f, IScriptEnvironment *env);
    void ready_data(Frame *f);
    void free_frame(Frame *f);
    void free_data(WorkingData *d);
};


class Filter : public GenericVideoFilter {
    PClip extmask;

    Parms Params;
    Buffer buf;
    FILE *outfile_ptr = nullptr;

    void find_segments(Frame *p, Frame *c, Frame *n);
    void find_motion_prev_cur(Frame *p, Frame *c);
    void find_motion_prev_next(Frame *p, Frame *c, Frame *n);

    void print_segments(int fn, IScriptEnvironment *env);

public:
    Filter(PClip _child, PClip _extmask, Parms _Params, IScriptEnvironment *env);
    ~Filter();
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment *env);
};


Filter::Filter(PClip _child, PClip _extmask, Parms _Params, IScriptEnvironment *env) :
    GenericVideoFilter(_child), extmask(_extmask), Params(_Params) {

    if (!vi.IsYV12() && !vi.IsYV16() && !vi.IsYV24())
        env->ThrowError("DeSpot: input to filter must be in YV12, YV16 or YV24!");

    buf.childclip = child;
    buf.num_frames = vi.num_frames;
    Params.height = vi.height;
    Params.width = vi.width;
    Params.pitch = ((vi.width + 7) / 8) * 8;
    Params.size = Params.height * Params.pitch;

    buf.data_cache[0].init(Params.size);
    buf.data_cache[1].init(Params.size);
    buf.motion = (BYTE *)malloc(Params.size);

    if (!Params.outfilename.empty()) {
        const double diag = sqrt(double(vi.width) * vi.width + double(vi.height) * vi.height);
        const int    thickness = int(floor(diag / 860 + 0.75));
        outfile_ptr = fopen(Params.outfilename.c_str(), "w");
        fprintf(outfile_ptr,
            "[Script Info]\n"
            "Title: DeSpot automatically generated file\n"
            "ScriptType: v4.00+\n"
            "WrapStyle: 0\n"
            "PlayResX: %d\n"
            "PlayResY: %d\n"
            "ScaledBorderAndShadow: yes\n"
            "Video Aspect Ratio: 0\n"
            "Video Zoom: 8\n"
            "Video Position: 0\n"
            "\n"
            "[V4+ Styles]\n"
            "Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding\n"
            "Style: Mask,Arial,20,&H00FFFFFF,&H000000FF,&H00000000,&H00000000,0,0,0,0,100,100,0,0,1,0,0,7,0,0,0,1\n"
            "Style: Outline,Arial,20,&HFFFFFFFF,&H000000FF,&H00FF00FF,&H00000000,0,0,0,0,100,100,0,0,1,%d,0,7,0,0,0,1\n"
            "\n"
            "[Events]\n"
            "Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text\n"
            "Comment: 0,0:00:00.00,0:00:00.00,Mask,,0,0,0,,--- Templates for manual editing ---\n"
            "Dialogue: 0,0:00:00.00,0:00:00.00,Mask,,0,0,0,,{\\pos(100,100)\\fscx100\\fscy100\\p1}m -4 -4 l 4 -4 4 4 -4 4{\\p0}\n"
            "Dialogue: 0,0:00:00.00,0:00:00.00,Mask,,0,0,0,,{\\pos(100,100)\\fscx100\\fscy100\\p1}m -6 -6 l 6 -6 6 6 -6 6{\\p0}\n"
            "Dialogue: 0,0:00:00.00,0:00:00.00,Mask,,0,0,0,,{\\pos(100,100)\\fscx100\\fscy100\\p1}m -8 -8 l 8 -8 8 8 -8 8{\\p0}\n"
            "Comment: 0,0:00:00.00,0:00:00.00,Mask,,0,0,0,,--- Masks ---\n",
            vi.width,
            vi.height,
            thickness
        );
    }
}

Filter::~Filter() {
    free(buf.motion);
    if (outfile_ptr != 0) {
        fclose(outfile_ptr);
        outfile_ptr = 0;
    }
}

/////////////////////////////////////////////////////////////////////
//
// Filter
//

void Filter::find_segments(Frame *p, Frame *c, Frame *n) {
    buf.alloc_noise(c);
    buf.alloc_segments(c);
    if (p && n) {
        ::find_outliers(p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, c->noise, Params);
        ::find_sizes(c->noise, c->segments, Params);
    }
}

void Filter::find_motion_prev_cur(Frame *p, Frame *c) {
    buf.alloc_motion(c);
    if (p && c) {
        ::find_motion(p->y, p->pitch, c->y, c->pitch, c->motion, Params);
    }

}

void Filter::find_motion_prev_next(Frame *p, Frame *c, Frame *n) {
    buf.alloc_motion(c);
    ::find_motion(p->y, p->pitch, n->y, n->pitch, c->motion, Params);
}

void	Filter::print_segments(int fn, IScriptEnvironment *env) {
    const int r = (Params.mc_flag) ? 3 : 1;
    const int fnd = fn / r;
    if (fn == fnd * r + ((r - 1) >> 1)) {
        Frame *c = buf.get_frame(fn, env); // fixme, what is this?
        const Segments &segs = c->segments;

        const int w = c->frame->GetRowSize();
        const int h = c->frame->GetHeight();

        const double fps = double(vi.fps_numerator) / double(vi.fps_denominator * r);

        ::print_segments(segs, outfile_ptr, fnd, fps, w, h, Params.spotmax1, Params.spotmax2);
    }
}

// ********************************************************************************

PVideoFrame __stdcall Filter::GetFrame(int fn, IScriptEnvironment *env) {
    Frame *pp = buf.get_frame(fn - 2, env);
    Frame *p = buf.get_frame(fn - 1, env);
    Frame *c = buf.get_frame(fn, env);
    Frame *n = buf.get_frame(fn + 1, env);
    Frame *nn = buf.get_frame(fn + 2, env);

    constexpr Exec *exec_median[3] = { cond_median, map_outliers, map_outliers };
    constexpr Exec *exec_pixel[3] = { remove_outliers, mark_outliers, map_outliers };

    if (p && n) {
        PVideoFrame fout = env->NewVideoFrame(vi);
        BYTE *fout_y = fout->GetWritePtr();
        int opitch = fout->GetPitch();

        if (Params.median) { // simple median mode (old)
            if (!c->motion) {
                find_motion_prev_cur(p, c);
                motion_denoise(c->motion, Params);
            }
            if (!n->motion) {
                find_motion_prev_cur(c, n);
                motion_denoise(n->motion, Params);
            }
            motion_merge(c->motion, n->motion, buf.motion, Params);
            (*exec_median[Params.show])(p->y, p->pitch, c->y, c->pitch, c->noise, buf.motion, n->y, n->pitch, fout_y, opitch, Params);
        } else if (Params.motpn && Params.seg != 0) { // motion prev-next,  segments mode
            // new mode - find motion prev to next
            if (!c->motion) {
                find_motion_prev_next(p, c, n);
                motion_denoise(c->motion, Params);
                motion_scene(c->motion, Params);
                if (extmask) {
                    PVideoFrame fextmask = extmask->GetFrame(fn, env);
                    add_external_mask(fextmask->GetReadPtr(), fextmask->GetPitch(), c->motion, Params.pitch, Params.width, Params.height);
                }
            }

            if (!c->noise) {
                find_segments(p, c, n);
                reject_on_motion(c->segments, c->motion, Params);
            }

            if (Params.show == 0) {
                vsh::bitblt(fout_y, opitch, c->y, c->pitch, Params.width, Params.height);
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                if (Params.tsmooth > 0)temporal_smooth(c->segments, p->y, p->pitch, n->y, n->pitch, fout_y, opitch, c->noise, c->motion, Params);
            } else if (Params.show == 1) {
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                mark_outliers(p->y, p->pitch, c->y, c->pitch, c->noise, c->motion, n->y, n->pitch, fout_y, opitch, Params);
            } else { // =2
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                map_outliers(p->y, p->pitch, c->y, c->pitch, c->noise, c->motion, n->y, n->pitch, fout_y, opitch, Params);
            }

            if (outfile_ptr != 0) {
                print_segments(fn, env);
            }
        } else if (Params.motpn && Params.seg == 0) { // motion prev-next, no segments mode
            // new mode - find motion prev to next
            if (!c->motion) {
                find_motion_prev_next(p, c, n);
                motion_denoise(c->motion, Params);
                motion_scene(c->motion, Params);
                if (extmask) {
                    PVideoFrame fextmask = extmask->GetFrame(fn, env);
                    add_external_mask(fextmask->GetReadPtr(), fextmask->GetPitch(), c->motion, Params.pitch, Params.width, Params.height);
                }
            }

            if (!c->noise) {
                find_segments(p, c, n);
                mark_noise(c->segments, c->noise, Params);
                noise_dilate(c->noise, Params);
            }


            (*exec_pixel[Params.show])(p->y, p->pitch, c->y, c->pitch, c->noise, c->motion, n->y, n->pitch, fout_y, opitch, Params);
        } else if (Params.seg != 0) { // motion prev-cur and cur-next (old mode), but segments mode
            if (!c->motion) {
                find_motion_prev_cur(p, c);
                find_motion_prev_cur(pp, p);
                if (!p->noise) {
                    find_segments(pp, p, c);
                    mark_noise(p->segments, p->noise, Params);
                }
                if (!c->noise) {
                    find_segments(p, c, n);
                    mark_noise(c->segments, c->noise, Params);
                }
                noise_to_one(p->noise, c->noise, c->motion, Params);
                motion_denoise(c->motion, Params);
                motion_scene(c->motion, Params);
            }

            if (!n->motion) {
                find_motion_prev_cur(c, n);
                if (!c->noise) {
                    find_segments(p, c, n);
                    mark_noise(c->segments, c->noise, Params);
                }
                if (!n->noise) {
                    find_segments(c, n, nn);
                    mark_noise(n->segments, n->noise, Params);
                }
                noise_to_one(c->noise, n->noise, n->motion, Params);
                motion_denoise(n->motion, Params);
                motion_scene(n->motion, Params);
            }
            motion_merge(c->motion, n->motion, buf.motion, Params);
            if (extmask) {
                PVideoFrame fextmask = extmask->GetFrame(fn, env);
                add_external_mask(fextmask->GetReadPtr(), fextmask->GetPitch(), buf.motion, Params.pitch, Params.width, Params.height);
            }
            reject_on_motion(c->segments, buf.motion, Params);
            if (Params.show == 0) {
                vsh::bitblt(fout_y, opitch, c->y, c->pitch, Params.width, Params.height); // copy as base
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                if (Params.tsmooth > 0) temporal_smooth(c->segments, p->y, p->pitch, n->y, n->pitch, fout_y, opitch, c->noise, c->motion, Params);
            } else if (Params.show == 1) {
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                mark_outliers(p->y, p->pitch, c->y, c->pitch, c->noise, buf.motion, n->y, n->pitch, fout_y, opitch, Params);
            } else { // =2
                remove_segments(c->segments, p->y, p->pitch, c->y, c->pitch, n->y, n->pitch, fout_y, opitch, c->noise, Params);
                map_outliers(p->y, p->pitch, c->y, c->pitch, c->noise, buf.motion, n->y, n->pitch, fout_y, opitch, Params);
            }
        } else if (Params.seg == 0) { // motion prev-cur and cur-next, no segments
            if (!c->motion) {
                find_motion_prev_cur(p, c);
                find_motion_prev_cur(pp, p);
                if (!p->noise) {
                    find_segments(pp, p, c);
                    mark_noise(p->segments, p->noise, Params);
                    noise_dilate(p->noise, Params);
                }
                if (!c->noise) {
                    find_segments(p, c, n);
                    mark_noise(c->segments, c->noise, Params);
                    noise_dilate(c->noise, Params);
                }
                noise_to_one(p->noise, c->noise, c->motion, Params);
                motion_denoise(c->motion, Params);
                motion_scene(c->motion, Params);
            }

            if (!n->motion) {
                find_motion_prev_cur(c, n);
                if (!c->noise) {
                    find_segments(p, c, n);
                    mark_noise(c->segments, c->noise, Params);
                    noise_dilate(c->noise, Params);
                }
                if (!n->noise) {
                    find_segments(c, n, nn);
                    mark_noise(n->segments, n->noise, Params);
                    noise_dilate(n->noise, Params);
                }
                noise_to_one(c->noise, n->noise, n->motion, Params);
                motion_denoise(n->motion, Params);
                motion_scene(n->motion, Params);
            }
            motion_merge(c->motion, n->motion, buf.motion, Params);
            if (extmask) {
                PVideoFrame fextmask = extmask->GetFrame(fn, env);
                add_external_mask(fextmask->GetReadPtr(), fextmask->GetPitch(), buf.motion, Params.pitch, Params.width, Params.height);
            }
            (*exec_pixel[Params.show])(p->y, p->pitch, c->y, c->pitch, c->noise, buf.motion, n->y, n->pitch, fout_y, opitch, Params);
        }

        // now process color planes
        if (Params.show == S_MAP && !Params.show_chroma) { //  grey map
            memset(fout->GetWritePtr(PLANAR_U), 0x80,
                fout->GetPitch(PLANAR_U) * fout->GetHeight(PLANAR_U));
            memset(fout->GetWritePtr(PLANAR_V), 0x80,
                fout->GetPitch(PLANAR_V) * fout->GetHeight(PLANAR_V));
        } else if (Params.show == S_MAP && Params.show_chroma) {// copy color planes
            vsh::bitblt(fout->GetWritePtr(PLANAR_U), fout->GetPitch(PLANAR_U),
                c->frame->GetReadPtr(PLANAR_U), c->frame->GetPitch(PLANAR_U),
                c->frame->GetRowSize(PLANAR_U), c->frame->GetHeight(PLANAR_U));
            vsh::bitblt(fout->GetWritePtr(PLANAR_V), fout->GetPitch(PLANAR_V),
                c->frame->GetReadPtr(PLANAR_V), c->frame->GetPitch(PLANAR_V),
                c->frame->GetRowSize(PLANAR_V), c->frame->GetHeight(PLANAR_V));
        } else if ((Params.show && S_MARK) && !Params.median) {// copy color planes
            vsh::bitblt(fout->GetWritePtr(PLANAR_U), fout->GetPitch(PLANAR_U),
                c->frame->GetReadPtr(PLANAR_U), c->frame->GetPitch(PLANAR_U),
                c->frame->GetRowSize(PLANAR_U), c->frame->GetHeight(PLANAR_U));
            vsh::bitblt(fout->GetWritePtr(PLANAR_V), fout->GetPitch(PLANAR_V),
                c->frame->GetReadPtr(PLANAR_V), c->frame->GetPitch(PLANAR_V),
                c->frame->GetRowSize(PLANAR_V), c->frame->GetHeight(PLANAR_V));
            // change color of marked noise to pink
            mark_color_plane(c->noise, fout->GetWritePtr(PLANAR_V), fout->GetPitch(PLANAR_V),
                fout->GetRowSize(PLANAR_V), fout->GetHeight(PLANAR_V), Params);
        } else if (Params.color && !Params.median) {	// normal mode with color
            // clean color YUV12 planes at places of luma spots (noise)
            clean_color_plane(p->frame->GetReadPtr(PLANAR_U), p->frame->GetPitch(PLANAR_U),
                c->frame->GetReadPtr(PLANAR_U), c->frame->GetPitch(PLANAR_U),
                n->frame->GetReadPtr(PLANAR_U), n->frame->GetPitch(PLANAR_U),
                c->frame->GetRowSize(PLANAR_U), c->frame->GetHeight(PLANAR_U),
                c->noise, fout->GetWritePtr(PLANAR_U), fout->GetPitch(PLANAR_U), Params);
            clean_color_plane(p->frame->GetReadPtr(PLANAR_V), p->frame->GetPitch(PLANAR_V),
                c->frame->GetReadPtr(PLANAR_V), c->frame->GetPitch(PLANAR_V),
                n->frame->GetReadPtr(PLANAR_V), n->frame->GetPitch(PLANAR_V),
                c->frame->GetRowSize(PLANAR_V), c->frame->GetHeight(PLANAR_V),
                c->noise, fout->GetWritePtr(PLANAR_V), fout->GetPitch(PLANAR_V), Params);
        } else { // mormal mode, copy color planes
            vsh::bitblt(fout->GetWritePtr(PLANAR_U), fout->GetPitch(PLANAR_U),
                c->frame->GetReadPtr(PLANAR_U), c->frame->GetPitch(PLANAR_U),
                c->frame->GetRowSize(PLANAR_U), c->frame->GetHeight(PLANAR_U));
            vsh::bitblt(fout->GetWritePtr(PLANAR_V), fout->GetPitch(PLANAR_V),
                c->frame->GetReadPtr(PLANAR_V), c->frame->GetPitch(PLANAR_V),
                c->frame->GetRowSize(PLANAR_V), c->frame->GetHeight(PLANAR_V));
        }

        return fout; // result
    } else { // first frame
        return c->frame; // nothing to do - no temporal info
    }
}


// macro for parameters defining
#define SET_INT(what, Min, Max) do {\
    i++;\
	tmp = args[i];\
    if (tmp.Defined()) {\
      int v = tmp.AsInt();\
      if (v < Min || v > Max)\
        env->ThrowError("DeSpot: "#what " must be from " #Min " to " #Max );\
      p.what = v;}} while (0)

AVSValue __cdecl Create_Despot(AVSValue args, void *user_data, IScriptEnvironment *env) {
    Parms p;
    AVSValue tmp;
    int i = 0;
    // default parameters values are set in despot.hpp
    SET_INT(mthres, 0, 255);
    SET_INT(mwidth, 0, 32000);
    SET_INT(mheight, 0, 32000);
    SET_INT(merode, 0, 100);
    i++;
    if (args[i].Defined()) p.y_next = args[i].AsBool(false) ? 2 : 1;  //interlaced
    i++;
    p.median = args[i].AsBool(false);
    SET_INT(p1, 0, 255);
    SET_INT(p2, 0, 255);
    SET_INT(pwidth, 0, 32000);
    SET_INT(pheight, 0, 32000);
    i++;
    p.ranked = args[i].AsBool(true);
    SET_INT(sign, -2, 2);
    SET_INT(maxpts, 0, 10000000);
    SET_INT(p1percent, 0, 100);
    SET_INT(dilate, 0, 1000);
    i++;
    p.fitluma = args[i].AsBool(false);
    SET_INT(blur, 0, 4);
    SET_INT(tsmooth, 0, 255);
    SET_INT(show, 0, 2);
    SET_INT(mark_v, 0, 255);
    i++;
    p.show_chroma = args[i].AsBool(false);
    i++;
    p.motpn = args[i].AsBool(true);
    SET_INT(seg, 0, 2);
    i++;
    p.color = args[i].AsBool(false);
    SET_INT(mscene, 0, 100);
    SET_INT(minpts, 0, 10000000);
    i++;
    int iextmask = i;
    i++; // planar
    p.outfilename = args[i].Defined() ? args[i].AsString() : "";
    ++i;
    p.mc_flag = args[i].AsBool(false);
    SET_INT(spotmax1, 1, 10000000);
    SET_INT(spotmax2, 1, 10000000);

    return new Filter(args[0].AsClip(), args[iextmask].Defined() ? args[iextmask].AsClip() : 0, p, env);
}

#define common_parms  "c[mthres]i[mwidth]i[mheight]i[merode]i[interlaced]b[median]b"
#define denoise_parms "[p1]i[p2]i[pwidth]i[pheight]i[ranked]b[sign]i[maxpts]i[p1percent]i[dilate]i[fitluma]b[blur]i[tsmooth]i[show]i[mark_v]i[show_chroma]b[motpn]b[seg]i[color]b[mscene]i[minpts]i[extmask]c[outfile]s[mc]b[spotmax1]i[spotmax2]i"

const AVS_Linkage *AVS_linkage = 0;
extern "C" __declspec(dllexport) const char *__stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage *const vectors) {
    AVS_linkage = vectors;
    env->AddFunction("DeSpot", common_parms denoise_parms, Create_Despot, 0);
    return "DeSpot plugin";
}

/////////////////////////////////////////////////////////////////////
//
// Frame
//

void WorkingData::init(size_t s) {
    noise = (BYTE *)malloc(s);
    motion = (BYTE *)malloc(s);
    segments = (Segment *)malloc(sizeof(Segment) * (s / 2 + 1));
}

WorkingData::~WorkingData() {
    if (noise)   free(noise);
    if (motion)  free(motion);
    if (segments) free(segments);
}

void Buffer::ready_frame(int n, Frame *f, IScriptEnvironment *env) {
    free_frame(f);
    f->num = n;
    f->frame = childclip->GetFrame(n, env);
    f->y = f->frame->GetReadPtr();
    f->pitch = f->frame->GetPitch();
}

void Buffer::ready_data(Frame *f) {
    WorkingData *d = data_cache + f->num % 2;
    free_data(d);
    d->num = f->num;
    f->data = d;
}

void Buffer::free_frame(Frame *f) {
    if (f->data) free_data(f->data);
    f->num = -1;
    f->frame = 0;
    f->y = 0;
}

void Buffer::free_data(WorkingData *d) {
    Frame *f = frame_cache + d->num % 16;
    if (f->num == d->num) {
        f->data = 0;
        f->noise = 0;
        f->motion = 0;
        f->segments.data = 0;
        f->segments.clear();
    }
    d->num = -1;
}

