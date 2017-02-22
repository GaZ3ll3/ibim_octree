//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_MOLECULE_H
#define IBIM_OCTREE_MOLECULE_H

#include "point.h"

class Molecule {
public:
    vector<point> centers;
    vector<scalar_t > radii;
    vector<scalar_t > charges;
    vector<vector<index_t >> adjacent;

    point center;
    scalar_t radius;
    index_t N;

    void load(std::string fileName) {
        std::ifstream pdbFile;
        std::string line;

        pdbFile.open(fileName);

        if (pdbFile.is_open()) {
            // read each line
            while (getline(pdbFile, line)) {
                if (line.find("ATOM") == 0) {
                    std::stringstream ss(line);
                    std::string buf;
                    vector<std::string> tokens;
                    while (ss >> buf) {
                        tokens.push_back(buf);
                    }

                    centers.push_back(point(
                            std::stod(tokens[5]),
                            std::stod(tokens[6]),
                            std::stod(tokens[7])
                    ));

                    radii.push_back(std::stof(tokens[9]));
                    charges.push_back(std::stof(tokens[8]));
                }
            }
            pdbFile.close();
            N = (index_t) centers.size();
            getCenter();
        }
    }
    void getAdjacent(scalar_t probe) {
        adjacent.resize(N);
        for (int i = 0; i < N; ++i) {
            scalar_t ri = radii[i] + probe;
            for (int j = 0; j < i; ++j) {
                // check if intersected.
                scalar_t  d = norm(centers[i] - centers[j]);
                scalar_t  rj = radii[j] + probe;
                if (d < ri + rj) {
                    adjacent[i].push_back(j);
                    adjacent[j].push_back(i);

                }
            }
        }
    }
    void getCenter() {
        assert(N > 0);
        point minP = {
                centers[0].data[0] + radii[0],
                centers[0].data[1] + radii[0],
                centers[0].data[2] + radii[0]};
        point maxP = {
                centers[0].data[0] - radii[0],
                centers[0].data[1] - radii[0],
                centers[0].data[2] - radii[0]};

        for (int i = 1; i < N; ++i) {
            point current_max_P = {
                    centers[i].data[0] + radii[i],
                    centers[i].data[1] + radii[i],
                    centers[i].data[2] + radii[i]
            };

            point current_min_P = {
                    centers[i].data[0] - radii[i],
                    centers[i].data[1] - radii[i],
                    centers[i].data[2] - radii[i]
            };


            for (int j = 0; j < 3; ++j) {
                if (maxP.data[j] <= current_max_P.data[j]) {maxP.data[j] = current_max_P.data[j];}

                if (minP.data[j] >= current_min_P.data[j]) {minP.data[j] = current_min_P.data[j];}
            }

        }

        this->center = (minP + maxP) * 0.5;
        this->radius = std::max(
                0.5 * (maxP.data[0] - minP.data[0]),
                std::max(
                        0.5 * (maxP.data[1] - minP.data[1]),
                        0.5 * (maxP.data[2] - minP.data[2]))
        );

    }


    scalar_t centralize(float range) {
        scalar_t s = range / radius;

        for (index_t i = 0; i < centers.size(); ++i) {
            centers[i].data[0] = s * (centers[i].data[0] - center.data[0]);
            centers[i].data[1] = s * (centers[i].data[1] - center.data[1]);
            centers[i].data[2] = s * (centers[i].data[2] - center.data[2]);
            radii[i] = s * radii[i];
        }
        return s;
    }
};


/*
 * float centerMolecule(Molecule mol, Molecule caTrace)
{
	// scale and translate the molecule so that the atom coordinates are between -250 and 250
	float minx,miny,minz,maxx,maxy,maxz;
	float s,sx,sy,sz;
	int i;

	minx = mol->xpoints[0];
	maxx = mol->xpoints[0];
	miny = mol->ypoints[0];
	maxy = mol->ypoints[0];
	minz = mol->zpoints[0];
	maxz = mol->zpoints[0];

	for (i=1;i<mol->npoints;i++){
		if (mol->xpoints[i]<minx)	minx=mol->xpoints[i];
		if (mol->xpoints[i]>maxx)	maxx=mol->xpoints[i];
		if (mol->ypoints[i]<miny)	miny=mol->ypoints[i];
		if (mol->ypoints[i]>maxy)	maxy=mol->ypoints[i];
		if (mol->zpoints[i]<minz)	minz=mol->zpoints[i];
		if (mol->zpoints[i]>maxz)	maxz=mol->zpoints[i];
	}
	sx = 400.0f/(maxx-minx);
	sy = 400.0f/(maxy-miny);
	sz = 400.0f/(maxz-minz);
	s = 0.0f;
	if (sx<sy && sx<sz)	s=sx;
	else if (sy<sx && sy<sz)s=sy;
	else s=sz;
	//printf("%f %f %f\n",sx,sy,sz);
	//printf("s = %f\n",s);
	printf("SPDBV quality = %f\n",(PR*s)/(512/N));

	for (i=0;i<mol->npoints;i++){
		mol->xpoints[i]=(s*(mol->xpoints[i]-((minx+maxx)/2)));
		mol->ypoints[i]=(s*(mol->ypoints[i]-((miny+maxy)/2)));
		mol->zpoints[i]=(s*(mol->zpoints[i]-((minz+maxz)/2)));
		mol->rpoints[i]=mol->rpoints[i]*s;
	}
	for (i=0;i<catrace->npoints;i++){
		catrace->xpoints[i]=(s*(caTrace->xpoints[i]-((minx+maxx)/2)));
		catrace->ypoints[i]=(s*(caTrace->ypoints[i]-((miny+maxy)/2)));
		catrace->zpoints[i]=(s*(caTrace->zpoints[i]-((minz+maxz)/2)));
		catrace->rpoints[i]=caTrace->rpoints[i]*s;
	}
	return s;
}
 */




#endif //IBIM_OCTREE_MOLECULE_H
