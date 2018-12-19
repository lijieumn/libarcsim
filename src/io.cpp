/*
  Copyright ©2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "io.hpp"

#include "util.hpp"
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <cassert>
#include <cfloat>
#include <json/json.h>
#include <fstream>
#include <png.h>
#include <sstream>
using namespace std;

// OBJ meshes

namespace arcsim {
    void get_valid_line(istream &in, string &line) {
        do
            getline(in, line);
        while (in && (line.length() == 0 || line[0] == '#'));
    }

    void triangle_to_obj(const string &inname, const string &outname) {
        std::cout << "outname is " << outname << std::endl;
        fstream outfile(outname.c_str(), ios::out);
        { // nodes
            string filename = inname + ".node";
            fstream file(filename.c_str(), ios::in);
            string line;
            get_valid_line(file, line);
            stringstream linestream(line);
            int nv, dim, na, nb;
            linestream >> nv >> dim >> na >> nb;
            for (int i = 0; i < nv; i++) {
                get_valid_line(file, line);
                stringstream linestream(line);
                int index;
                linestream >> index;
                Vec2 u;
                linestream >> u[0] >> u[1];
                outfile << "v " << u[0] << " " << u[1] << " " << 0 << endl;
            }
        }
        { // eles
            string filename = inname + ".ele";
            fstream file(filename.c_str(), ios::in);
            string line;
            get_valid_line(file, line);
            stringstream linestream(line);
            int nt, nn, na;
            linestream >> nt >> nn >> na;
            for (int i = 0; i < nt; i++) {
                get_valid_line(file, line);
                stringstream linestream(line);
                int index;
                linestream >> index;
                int v0, v1, v2;
                linestream >> v0 >> v1 >> v2;
                outfile << "f " << v0 + 1 << " " << v1 + 1 << " " << v2 + 1 << endl;
            }
        }
    }


    void load_vdb(const string &filename) {

        VisualDebugger *vd = VisualDebugger::getInstance();

        fstream file(filename.c_str(), ios::in);
        if (!file) {
            cout << "Error: failed to open file " << filename << endl;
            return;
        }

        Vec3 color;
        unsigned char key;

        while (file) {
            string line;
            get_valid_line(file, line);
            stringstream linestream(line);
            string keyword;
            linestream >> keyword;
            if (keyword == "vis") {
                linestream >> key >> color[0] >> color[1] >> color[2];
                //vd->infoList.push_back( VisualInfo(key,color) );
                //vd->keys.push_back( key );
                //vd->display.push_back( false );
            } else if (keyword == "pos") {
                Vec3 pos;
                linestream >> pos[0] >> pos[1] >> pos[2];
                vd->addVisualPoint3(pos, color, key);
            } else if (keyword == "upos") {
                Vec2 upos;
                linestream >> upos[0] >> upos[1];
                vd->addVisualPoint2(upos, color, key);
            } else if (keyword == "line") {
                Vec3 pos0;
                Vec3 pos1;
                linestream >> pos0[0] >> pos0[1] >> pos0[2] >> pos1[0] >> pos1[1] >> pos1[2];
                vd->addVisualLine3(pos0, pos1, color, key);
            } else if (keyword == "uline") {
                Vec2 upos0;
                Vec2 upos1;
                linestream >> upos0[0] >> upos0[1] >> upos1[0] >> upos1[1];
                vd->addVisualLine2(upos0, upos1, color, key);
            }
        }


    }

    void save_vdb(const std::string &filename) {


        VisualDebugger *vd = VisualDebugger::getInstance();
        fstream file(filename.c_str(), ios::out);
        file.precision(17);

        for (int i = 0; i < vd->infoList.size(); i++) {
            VisualInfo &visInfo = vd->infoList[i];
            file << "vis" << " " << visInfo.key << " " << visInfo.color[0] << " " << visInfo.color[1] << " "
                 << visInfo.color[2] << std::endl;
            for (int j = 0; j < visInfo.pos.size(); j++) {
                file << "pos" << " " << visInfo.pos[j][0] << " " << visInfo.pos[j][1] << " " << visInfo.pos[j][2]
                     << std::endl;
            }

            for (int j = 0; j < visInfo.u_pos.size(); j++) {
                file << "upos" << " " << visInfo.u_pos[j][0] << " " << visInfo.u_pos[j][1] << " " << visInfo.u_pos[j][2]
                     << std::endl;
            }

            for (int j = 0; j < visInfo.lines.size(); j++) {
                Vec3 firstPt = visInfo.lines[j].first;
                Vec3 secondPt = visInfo.lines[j].second;
                file << "line" << " " << firstPt[0] << " " << firstPt[1] << " " << firstPt[2] << " "
                     << secondPt[0] << " " << secondPt[1] << " " << secondPt[2] << std::endl;
            }

            for (int j = 0; j < visInfo.u_lines.size(); j++) {
                Vec2 firstPt = visInfo.u_lines[j].first;
                Vec2 secondPt = visInfo.u_lines[j].second;
                file << "uline" << " " << firstPt[0] << " " << firstPt[1] << " " << secondPt[0] << " " << secondPt[1]
                     << std::endl;
            }

        }

        file.close();

    }


    void load_obj(Mesh &mesh, const string &filename) {

        delete_mesh(mesh);
        fstream file(filename.c_str(), ios::in);
        if (!file) {
            cout << "Error: failed to open file " << filename << endl;
            return;
        }
        bool load_adjf = false;
        bool load_adje = false;
        double scale = 1;
        while (file) {
            string line;
            get_valid_line(file, line);
            stringstream linestream(line);
            string keyword;
            linestream >> keyword;
            if (keyword == "scale") {
                linestream >> scale;
            } else if (keyword == "vt") {
                Vec2 u;
                linestream >> u[0] >> u[1];
                u[0] *= scale;
                u[1] *= scale;
                mesh.add(new Vert(u, -1));
            } else if (keyword == "vl") {
                linestream >> mesh.verts.back()->label;
            } else if (keyword == "adjf") {
                load_adjf = true;
                vector<unsigned int> &indices = mesh.verts.back()->adjf_indices;
                unsigned int index;
                while (linestream >> index) {
                    indices.push_back(index);
                }
            } else if (keyword == "v") {
                Vec3 x;
                linestream >> x[0] >> x[1] >> x[2];
                mesh.add(new Node(x, Vec3(0)));
            } else if (keyword == "temp") {
                mesh.nodes.back()->temp = true;
            } else if (keyword == "temp2") {
                mesh.nodes.back()->temp2 = true;
            } else if (keyword == "vv") {
                Vec3 &v = mesh.nodes.back()->v;
                linestream >> v[0] >> v[1] >> v[2];
            } else if (keyword == "ny") {
                Vec3 &y = mesh.nodes.back()->y;
                linestream >> y[0] >> y[1] >> y[2];
            } else if (keyword == "nv") {
                Vec3 &v = mesh.nodes.back()->v;
                linestream >> v[0] >> v[1] >> v[2];
            } else if (keyword == "nl") {
                linestream >> mesh.nodes.back()->label;
            } else if (keyword == "nf") {
                Vec3 &r = mesh.nodes.back()->r;
                linestream >> r[0] >> r[1] >> r[2];
            } else if (keyword == "adje") {
                load_adje = true;
                vector<unsigned int> &indices = mesh.nodes.back()->adje_indices;
                unsigned int index;
                while (linestream >> index) {
                    indices.push_back(index);
                }
            } else if (keyword == "e") {
                int n0, n1;
                linestream >> n0 >> n1;
                mesh.add(new Edge(mesh.nodes[n0 - 1], mesh.nodes[n1 - 1]));
            } else if (keyword == "ea") {
                linestream >> mesh.edges.back()->theta_ideal;
            } else if (keyword == "ed") {
                linestream >> mesh.edges.back()->damage;
            } else if (keyword == "el") {
                linestream >> mesh.edges.back()->label;
            } else if (keyword == "f") {
                vector<Vert *> verts;
                vector<Node *> nodes;
                string w;
                while (linestream >> w) {
                    stringstream wstream(w);
                    int v, n;
                    char c;
                    wstream >> n >> c >> v;
                    nodes.push_back(mesh.nodes[n - 1]);
                    if (wstream) {
                        verts.push_back(mesh.verts[v - 1]);
                    } else if (!nodes.back()->verts.empty()) {
                        verts.push_back(nodes.back()->verts[0]);
                    } else {
                        verts.push_back(new Vert(project<2>(nodes.back()->x),
                                                 nodes.back()->label));
                        mesh.add(verts.back());
                    }
                }
                for (int v = 0; v < verts.size(); v++)
                    connect(verts[v], nodes[v]);
                vector<Face *> faces = triangulate(verts);
                for (int f = 0; f < faces.size(); f++)
                    mesh.add(faces[f]);
            } else if (keyword == "tl" || keyword == "fl") {
                linestream >> mesh.faces.back()->label;
            } else if (keyword == "ts" || keyword == "fs") {
                Mat2x2 &S = mesh.faces.back()->S_plastic;
                linestream >> S(0, 0) >> S(0, 1) >> S(1, 0) >> S(1, 1);
            } else if (keyword == "td" || keyword == "fd") {
                linestream >> mesh.faces.back()->damage;
            }
        }
        if (load_adjf) {
            for (int v = 0; v < mesh.verts.size(); v++) {
                Vert *vert = mesh.verts[v];
                vert->adjf.clear();
                for (int f = 0; f < vert->adjf_indices.size(); f++) {
                    unsigned int index = vert->adjf_indices[f];
                    vert->adjf.push_back(mesh.faces[index]);
                }
            }
        }
        if (load_adje) {
            for (int n = 0; n < mesh.nodes.size(); n++) {
                Node *node = mesh.nodes[n];
                node->adje.clear();
                for (int e = 0; e < node->adje_indices.size(); e++) {
                    unsigned int index = node->adje_indices[e];
                    node->adje.push_back(mesh.edges[index]);
                }
            }
        }
        mark_nodes_to_preserve(mesh);
        compute_ms_data(mesh);
    }

    void load_objs(vector<Mesh *> &meshes, const string &prefix) {
        for (int m = 0; m < meshes.size(); m++)
            load_obj(*meshes[m], stringf("%s_%02d.obj", prefix.c_str(), m));
    }

    void load_obs_frames(vector<Mesh> &meshes, const string &prefix, const int frameStart, const int frameEnd) {
        meshes.resize(frameEnd - frameStart + 1);
        for (int m = frameStart; m <= frameEnd; m++) {
            load_obj(meshes[m], stringf("%s_%04d.obj", prefix.c_str(), m));
        }
    }


    static double angle(const Vec3 &x0, const Vec3 &x1, const Vec3 &x2) {
        Vec3 e1 = normalize(x1 - x0);
        Vec3 e2 = normalize(x2 - x0);
        return acos(clamp(dot(e1, e2), -1., 1.));
    }

    vector<Face *> triangulate(const vector<Vert *> &verts) {
        int n = verts.size();
        double best_min_angle = 0;
        int best_root = -1;
        for (int i = 0; i < n; i++) {
            double min_angle = infinity;
            const Vert *vert0 = verts[i];
            for (int j = 2; j < n; j++) {
                const Vert *vert1 = verts[(i + j - 1) % n], *vert2 = verts[(i + j) % n];
                min_angle = min(min_angle,
                                angle(vert0->node->x, vert1->node->x, vert2->node->x),
                                angle(vert1->node->x, vert2->node->x, vert0->node->x),
                                angle(vert2->node->x, vert0->node->x, vert1->node->x));
            }
            // if (min_angle > best_min_angle) {
            best_min_angle = min_angle;
            best_root = i;
            // }
        }
        int i = best_root;
        Vert *vert0 = verts[i];
        vector<Face *> tris;
        for (int j = 2; j < n; j++) {
            Vert *vert1 = verts[(i + j - 1) % n], *vert2 = verts[(i + j) % n];
            tris.push_back(new Face(vert0, vert1, vert2, 0));
        }
        return tris;
    }

    void save_obj(const Mesh &mesh, const string &filename) {
        fstream file(filename.c_str(), ios::out);
        file.precision(17);
        Vec2 min, max;
        for (int v = 0; v < mesh.verts.size(); v++) {
            const Vert *vert = mesh.verts[v];
            if (v == 0) {
                min = vert->u;
                max = vert->u;
            } else {
                if (min[0] > vert->u[0]) {
                    min[0] = vert->u[0];
                }
                if (min[1] > vert->u[1]) {
                    min[1] = vert->u[1];
                }
                if (max[0] < vert->u[0]) {
                    max[0] = vert->u[0];
                }
                if (max[1] < vert->u[1]) {
                    max[1] = vert->u[1];
                }
            }
        }
        Vec2 size = max - min;
        double scale = size[0] > size[1] ? size[0] : size[1];
        Vec2 shift = min;
        double diff = fabs(size[0] - size[1]) / 2;
        if (size[0] > size[1]) {
            shift[1] -= diff;
        } else {
            shift[0] -= diff;
        }
        file << "scale " << scale << endl;
        for (int v = 0; v < mesh.verts.size(); v++) {
            const Vert *vert = mesh.verts[v];
            Vec2 u = vert->u - shift;
            file << "vt " << u[0] / scale << " " << u[1] / scale << endl;
            if (vert->label)
                file << "vl " << vert->label << endl;
            if (vert->adjf.size() > 0) {
                file << "adjf ";
                for (int a = 0; a < vert->adjf.size(); a++) {
                    file << vert->adjf[a]->index << " ";
                }
                file << endl;
            }
        }
        for (int n = 0; n < mesh.nodes.size(); n++) {
            const Node *node = mesh.nodes[n];
            file << "v " << node->x[0] << " " << node->x[1] << " "
                 << node->x[2] << endl;
            if (node->temp)
                file << "temp" << endl;
            if (node->temp2)
                file << "temp2" << endl;
            if (norm2(node->x - node->y))
                file << "ny " << node->y[0] << " " << node->y[1] << " "
                     << node->y[2] << endl;
            if (norm2(node->v))
                file << "nv " << node->v[0] << " " << node->v[1] << " "
                     << node->v[2] << endl;
            if (node->label)
                file << "nl " << node->label << endl;
            if (norm2(node->r))
                file << "nf " << node->r[0] << " " << node->r[1] << " "
                     << node->r[2] << endl;
            if (node->adje.size() > 0) {
                file << "adje ";
                for (int a = 0; a < node->adje.size(); a++) {
                    file << node->adje[a]->index << " ";
                }
                file << endl;
            }
        }
        for (int e = 0; e < mesh.edges.size(); e++) {
            const Edge *edge = mesh.edges[e];
            // if (edge->theta_ideal || edge->label) {
            file << "e " << edge->n[0]->index + 1 << " " << edge->n[1]->index + 1
                 << endl;
            if (edge->theta_ideal)
                file << "ea " << edge->theta_ideal << endl;
            if (edge->damage)
                file << "ed " << edge->damage << endl;
            if (edge->label)
                file << "el " << edge->label << endl;
            // }
        }
        for (int f = 0; f < mesh.faces.size(); f++) {
            const Face *face = mesh.faces[f];
            file << "f " << face->v[0]->node->index + 1 << "/" << face->v[0]->index + 1
                 << " " << face->v[1]->node->index + 1 << "/" << face->v[1]->index + 1
                 << " " << face->v[2]->node->index + 1 << "/" << face->v[2]->index + 1
                 << endl;
            if (face->label)
                file << "tl " << face->label << endl;
            if (norm2_F(face->S_plastic)) {
                const Mat2x2 &S = face->S_plastic;
                file << "ts " << S(0, 0) << " " << S(0, 1) << " " << S(1, 0) << " "
                     << S(1, 1) << endl;
            }
            if (face->damage)
                file << "td " << face->damage << endl;
        }
    }

    void save_objs(const vector<Mesh *> &meshes, const string &prefix) {
        for (int m = 0; m < meshes.size(); m++)
            save_obj(*meshes[m], stringf("%s_%02d.obj", prefix.c_str(), m));
    }

    void save_transformation(const Transformation &tr, const string &filename) {
        FILE *file = fopen(filename.c_str(), "w");
        pair<Vec3, double> axis_angle = tr.rotation.to_axisangle();
        Vec3 axis = axis_angle.first;
        double angle = axis_angle.second * 180 / M_PI;
        fprintf(file, "<rotate angle=\"%f\" x=\"%f\" y=\"%f\" z=\"%f\"/>\n",
                angle, axis[0], axis[1], axis[2]);
        fprintf(file, "<scale value=\"%f\"/>\n", tr.scale);
        fprintf(file, "<translate x=\"%f\" y=\"%f\" z=\"%f\"/>\n",
                tr.translation[0], tr.translation[1], tr.translation[2]);
        fclose(file);
    }

// Images

    void flip_image(int w, int h, unsigned char *pixels);

    void save_png(const char *filename, int width, int height,
                  unsigned char *pixels, bool has_alpha = false);

    void flip_image(int w, int h, unsigned char *pixels) {
        for (int j = 0; j < h / 2; j++)
            for (int i = 0; i < w; i++)
                for (int c = 0; c < 3; c++)
                    swap(pixels[(i + w * j) * 3 + c], pixels[(i + w * (h - 1 - j)) * 3 + c]);
    }

    void save_png(const char *filename, int width, int height,
                  unsigned char *pixels, bool has_alpha) {
#ifndef _WIN32
        FILE *file = fopen(filename, "wb");
        if (!file) {
            printf("Couldn't open file %s for writing.\n", filename);
            return;
        }
        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
                                                      NULL, NULL);
        if (!png_ptr) {
            printf("Couldn't create a PNG write structure.\n");
            fclose(file);
            return;
        }
        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            printf("Couldn't create a PNG info structure.\n");
            png_destroy_write_struct(&png_ptr, NULL);
            fclose(file);
            return;
        }
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Had a problem writing %s.\n", filename);
            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(file);
            return;
        }
        png_init_io(png_ptr, file);
        png_set_IHDR(png_ptr, info_ptr, width, height, 8,
                     has_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
                     PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                     PNG_FILTER_TYPE_DEFAULT);
        int channels = has_alpha ? 4 : 3;
        png_bytep *row_pointers = (png_bytep *) new unsigned char *[height];
        for (int y = 0; y < height; y++)
            row_pointers[y] = (png_bytep) &pixels[y * width * channels];
        png_set_rows(png_ptr, info_ptr, row_pointers);
        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
        delete[] row_pointers;
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(file);
#endif
    }

    void ensure_existing_directory(const std::string &path) {
        using namespace boost::filesystem;
        if (!exists(path))
            create_directory(path);
        if (!is_directory(path)) {
            cout << "Error: " << path << " is not a directory!" << endl;
            abort();
        }
    }


    void copy_file(const string &input, const string &output) {
        if (input == output) {
            return;
        }
        if (boost::filesystem::exists(output)) {
            boost::filesystem::remove(output);
        }
        boost::filesystem::copy_file(
                input, output);
    }
}