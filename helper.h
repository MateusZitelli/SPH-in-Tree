//SPH in Tree - SPH using barnes-hut algorithm

//Copyright 2012 Mateus Zitelli <zitellimateus@gmail.com>

//This program is free software; you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation; either version 2 of the License, or
//(at your option) any later version.

//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//MA 02110-1301, USA.

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "./Eigen/Dense"
#define MAX_NODES 150000
#define BODIES_QUANTITY 800
#define WIDTH 800
#define HEIGHT 800
#define PI 3.141592
#define SIZE_OF_SIMULATION 500
#define ALPHA 0.5
#define BETA 2 * ALPHA

int node_quantity;
int roots_quantity;
FILE *positionData;

using namespace Eigen;
using namespace std;

typedef Matrix<double, 3, 1> vector3d;

struct body {
	vector3d position;
	vector3d force;
	vector3d speed;
	double density;
	double pression;
	double h;
	double Ss;
	double acel;
	double mass;
};

struct node {
	vector3d start;
	vector3d end;
	int bodies_quantity;
	int initialized;
	int has_center_of_mass;
	int deep;
	struct body **bodies;
	struct node *UNE;
	struct node *USE;
	struct node *USW;
	struct node *UNW;
	struct node *DNE;
	struct node *DSE;
	struct node *DSW;
	struct node *DNW;
	struct node *UP;
	struct body centerOfMass;
};

struct node *nodes;
struct node **roots;
struct body *bodies;

double W(double r, double h){
        double R = r/h;
        if(R > 1) return(0);
        double fact = 8.0 / (PI * h * h * h);
        if(R >= 0 && R <= 0.5){
                double R2 = R * R;
                return(fact * (1 - 6 * R2 + 6 * R2 * R));
        }else if(R > 0.5 && R <= 1){
                return(fact * 2 * (1 - R) * (1 - R) * (1 - R));
        }
}

void initializeNode(struct node *node, struct node *up, double sx, double sy,double sz, double ex, double ey, double ez,int bodies_quantity, int deep){
	vector3d coord;
	float regulator;
	coord(0) = sx;
	coord(1) = sy;
	coord(2) = sz;
	node->start = coord;
	coord(0) = ex;
	coord(1) = ey;
	coord(2) = ez;
	node->end = coord;
	node->deep = deep;
	node->bodies_quantity = 0;
	node->has_center_of_mass = 0;
	if (!node->initialized) {
		node->bodies =
		    (struct body **)malloc(sizeof(struct body *) *
					   bodies_quantity);
		node->initialized = 1;
	}
	node->UNE = NULL;
	node->UNW = NULL;
	node->USE = NULL;
	node->USW = NULL;
	node->DNE = NULL;
	node->DNW = NULL;
	node->DSE = NULL;
	node->DSW = NULL;
	node->UP = up;
}

void addBodyInNode(struct body *body, struct node *node)
{
	node->bodies[node->bodies_quantity++] = body;
}

void resetNodes(void)
{
        nodes[0].UNE = NULL;
        nodes[0].UNW = NULL;
        nodes[0].USE = NULL;
        nodes[0].USW = NULL;
        nodes[0].DNE = NULL;
        nodes[0].DNW = NULL;
        nodes[0].DSE = NULL;
        nodes[0].DSW = NULL;
        nodes[0].UP = NULL;
        nodes[0].initialized = 0;
        nodes[0].has_center_of_mass = 0;
        int i;
        printf("Reset %i\n", node_quantity);
        for(i = 1;i < node_quantity; i++){
                nodes[i].UNE = NULL;
                nodes[i].UNW = NULL;
                nodes[i].USE = NULL;
                nodes[i].USW = NULL;
                nodes[i].DNE = NULL;
                nodes[i].DNW = NULL;
                nodes[i].DSE = NULL;
                nodes[i].DSW = NULL;
                nodes[i].UP = NULL;
                nodes[i].has_center_of_mass = 0;
                if(nodes[i].initialized && nodes[i].bodies != NULL){
	               free(nodes[i].bodies);
                }
	        nodes[i].initialized = 0;
        }
        node_quantity = 1;
}

void divideNode(struct node *node)
{
	if (node == NULL)
		return;
	if (node->bodies_quantity == 1 && node->bodies_quantity != 0) {
		roots[roots_quantity++] = node;
		return;
	} else if (node->bodies_quantity == 0) {
		return;
	}
	int i, UNW = 0, UNE = 0, USW = 0, USE = 0, DNW = 0, DNE = 0, DSW =
	    0, DSE = 0;
	double sx, sy, sz, ex, ey, ez, mx, my, mz;
	sx = node->start(0);
	sy = node->start(1);
	sz = node->start(2);
	ex = node->end(0);
	ey = node->end(1);
	ez = node->end(2);
	mx = (sx + ex) / 2.0;
	my = (sy + ey) / 2.0;
	mz = (sz + ez) / 2.0;
	for (i = 0; i < node->bodies_quantity; i++) {
		if (node->bodies[i]->position(0) < sx
		    || node->bodies[i]->position(0) >= ex
		    || node->bodies[i]->position(1) < sy
		    || node->bodies[i]->position(1) >= ey
		    || node->bodies[i]->position(2) < sz
		    || node->bodies[i]->position(2) >= ez)
			continue;
		if (node->bodies[i]->position(0) < mx) {
			if (node->bodies[i]->position(1) < my) {
				if (node->bodies[i]->position(2) < mz) {
					++UNW;
				} else {	// z >= mz
					++DNW;
				}
			} else {	// y >= my
				if (node->bodies[i]->position(2) < mz) {
					++USW;
				} else {	// z >= mz
					++DSW;
				}
			}
		} else {	// x >= mx
			if (node->bodies[i]->position(1) < my) {
				if (node->bodies[i]->position(2) < mz) {
					++UNE;
				} else {	// z >= mz
					++DNE;
				}
			} else {	// y >= my
				if (node->bodies[i]->position(2) < mz) {
					++USE;
				} else {	// z >= mz
					++DSE;
				}
			}
		}
	}
	for (i = 0; i < node->bodies_quantity; i++) {
		if (node->bodies[i]->position(0) < sx
		    || node->bodies[i]->position(0) >= ex
		    || node->bodies[i]->position(1) < sy
		    || node->bodies[i]->position(1) >= ey
		    || node->bodies[i]->position(2) < sz
		    || node->bodies[i]->position(2) >= ez)
			continue;
		if (node->bodies[i]->position(0) < mx) {
			if (node->bodies[i]->position(1) < my) {
				if (node->bodies[i]->position(2) < mz) {
					if (node->UNW == NULL) {
						node->UNW =
						    &nodes[node_quantity++];
						initializeNode(node->UNW, node,
							       sx, sy, sz, mx,
							       my, mz, UNW,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->UNW);
				} else {	// z >= mz
					if (node->DNW == NULL) {
						node->DNW =
						    &nodes[node_quantity++];
						initializeNode(node->DNW, node,
							       sx, sy, mz, mx,
							       my, ez, DNW,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->DNW);
				}
			} else {	// y >= my
				if (node->bodies[i]->position(2) < mz) {
					if (node->USW == NULL) {
						node->USW =
						    &nodes[node_quantity++];
						initializeNode(node->USW, node,
							       sx, my, sz, mx,
							       ey, mz, USW,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->USW);
				} else {	// z >= mz
					if (node->DSW == NULL) {
						node->DSW =
						    &nodes[node_quantity++];
						initializeNode(node->DSW, node,
							       sx, my, mz, mx,
							       ey, ez, DSW,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->DSW);
				}
			}
		} else {	// x >= mx
			if (node->bodies[i]->position(1) < my) {
				if (node->bodies[i]->position(2) < mz) {
					if (node->UNE == NULL) {
						node->UNE =
						    &nodes[node_quantity++];
						initializeNode(node->UNE, node,
							       mx, sy, sz, ex,
							       my, mz, UNE,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->UNE);
				} else {	// z >= mz
					if (node->DNE == NULL) {
						node->DNE =
						    &nodes[node_quantity++];
						initializeNode(node->DNE, node,
							       mx, sy, mz, ex,
							       my, ez, DNE,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->DNE);
				}
			} else {	// y >= my
				if (node->bodies[i]->position(2) < mz) {
					if (node->USE == NULL) {
						node->USE =
						    &nodes[node_quantity++];
						initializeNode(node->USE, node,
							       mx, my, sz, ex,
							       ey, mz, USE,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->USE);
				} else {	// z >= mz
					if (node->DSE == NULL) {
						node->DSE =
						    &nodes[node_quantity++];
						initializeNode(node->DSE, node,
							       mx, my, mz, ex,
							       ey, ez, DSE,
							       node->deep + 1);
					}
					addBodyInNode(node->bodies[i],
						      node->DSE);
				}
			}
		}
	}
	divideNode(node->UNW);
	divideNode(node->UNE);
	divideNode(node->USW);
	divideNode(node->USE);
	divideNode(node->DNW);
	divideNode(node->DNE);
	divideNode(node->DSW);
	divideNode(node->DSE);
}

void setCenterOfMass(struct node *node)
{
	double px = 0, py = 0, pz = 0, mass = 0, totalSpeed, relativisticAjust;
	double vx = 0, vy = 0, vz = 0, d= 0, h = 0;
	int i;
	struct body centerOfMass;
	if (node->bodies_quantity >= 1) {
		for (i = 0; i < node->bodies_quantity; i++) {
			px +=
			    node->bodies[i]->position(0) *
			    node->bodies[i]->mass;
			py +=
			    node->bodies[i]->position(1) *
			    node->bodies[i]->mass;
			pz +=
			    node->bodies[i]->position(2) *
			    node->bodies[i]->mass;
			vx +=
			    node->bodies[i]->speed(0) *
			    node->bodies[i]->mass;
			vy +=
			    node->bodies[i]->speed(1) *
			    node->bodies[i]->mass;
			vz +=
			    node->bodies[i]->speed(2) *
			    node->bodies[i]->mass;
			d += node->bodies[i]->density * node->bodies[i]->mass;
			mass += node->bodies[i]->mass;
			h += node->bodies[i]->h * node->bodies[i]->mass;
		}
		centerOfMass.position(0) = px / mass;
		centerOfMass.position(1) = py / mass;
		centerOfMass.position(2) = pz / mass;
		centerOfMass.speed(0) = vx / mass;
		centerOfMass.speed(1) = vy / mass;
		centerOfMass.speed(2) = vz / mass;
		centerOfMass.mass = mass;
		centerOfMass.density = d / mass;
		centerOfMass.h = sqrt(h / mass);
		node->centerOfMass = centerOfMass;
	} else {
		node->centerOfMass.position(0) = 0;
		node->centerOfMass.position(1) = 0;
		node->centerOfMass.position(2) = 0;
		node->centerOfMass.mass = 0;
		node->centerOfMass.density = 0;
		node->centerOfMass.h = 0;
	}
	node->has_center_of_mass = 1;
}
/*
void applyForceBetweenBodies(struct body *b1, struct body *b2)
{
	double xDistance, yDistance, zDistance, DistanceSquared, force,
	    distance;
	xDistance = b2->position(0) - b1->position(0);
	yDistance = b2->position(1) - b1->position(1);
	zDistance = b2->position(2) - b1->position(2);
	DistanceSquared =
	    xDistance * xDistance + yDistance * yDistance +
	    zDistance * zDistance + EPS2;
	force = K * b2->mass * b1->mass / DistanceSquared ;
	double dist = sqrt(DistanceSquared);
	b1->force(0) += xDistance / dist * force;
	b1->force(1) += yDistance / dist * force;
	b1->force(2) += zDistance / dist * force;
}
*/
int applyForce(struct node *node, struct body *body)
{
	double xDistance, yDistance, zDistance, DistanceSquared, force,
	    distance;
	struct body centerOfMass;
	xDistance = (body->position(0) - (node->end(0) + node->start(0)) / 2.0);
	yDistance = (body->position(1) - (node->end(1) + node->start(1)) / 2.0);
	zDistance = (body->position(2) - (node->end(2) + node->start(2)) / 2.0);
	distance =
	    sqrt(xDistance * xDistance + yDistance * yDistance +
		 zDistance * zDistance);
	if ((node->end(0) - node->start(0)) / distance < ALPHA
	    || node->bodies_quantity == 1) {
		centerOfMass = node->centerOfMass;
		vector3d diffSpeed = (body->speed - centerOfMass.speed);
                vector3d diffPosition(xDistance, yDistance, zDistance);
                double dotP = diffSpeed.dot(diffPosition);
                double visocity_factor = 0;
                if(dotP < 0){
                        double mi = ((body->h + centerOfMass.h) / 2.0) * dotP / (diffPosition.squaredNorm());
                        visocity_factor = (-ALPHA * ((body->Ss + centerOfMass.Ss) / 2.0) * mi + BETA * mi * mi) / ((body->density + centerOfMass.density) / 2.0);
                }
                double force = centerOfMass.mass * body->mass * (centerOfMass.pression / (centerOfMass.density * centerOfMass.density) + body->pression / (body->density * body->density)) * W(distance, centerOfMass.h);
                double vel = body->mass * visocity_factor * diffSpeed.norm() * W(distance , ((body->h + centerOfMass.h) / 2.0));
                body->force(0) += force * xDistance / distance;
                body->force(1) += force * yDistance / distance;
                body->force(2) += force * zDistance / distance;
                body->speed(0) +=  10 *vel * xDistance / distance;
                body->speed(1) +=  10 *vel * yDistance / distance;
                body->speed(2) +=  10 *vel * zDistance / distance;
                //cout << force << "\n";
		return (1);
	} else {
		return (0);
	}
}

void forceOverNode(struct node *node, struct node *down, struct body *body,
		   int inverse)
{
	int var, i;
	if (node == NULL)
		return;
	if (node->UNE != NULL && (node->UNE != down || inverse)) {
		var = 1;
		var = applyForce(node->UNE, body);
		if (!var) {
			forceOverNode(node->UNE, node, body, 1);
		}
	}
	if (node->UNW != NULL && (node->UNW != down || inverse)) {
		var = applyForce(node->UNW, body);
		if (!var) {
			forceOverNode(node->UNW, node, body, 1);
		}
	}
	if (node->USE != NULL && (node->USE != down || inverse)) {
		var = applyForce(node->USE, body);
		if (!var) {
			forceOverNode(node->USE, node, body, 1);
		}
	}
	if (node->USW != NULL && (node->USW != down || inverse)) {
		var = applyForce(node->USW, body);
		if (!var) {
			forceOverNode(node->USW, node, body, 1);
		}
	}
	if (node->DNE != NULL && (node->DNE != down || inverse)) {
		var = applyForce(node->DNE, body);
		if (!var) {
			forceOverNode(node->DNE, node, body, 1);
		}
	}
	if (node->DNW != NULL && (node->DNW != down || inverse)) {
		var = applyForce(node->DNW, body);
		if (!var) {
			forceOverNode(node->DNW, node, body, 1);
		}
	}
	if (node->DSE != NULL && (node->DSE != down || inverse)) {
		var = applyForce(node->DSE, body);
		if (!var) {
			forceOverNode(node->DSE, node, body, 1);
		}
	}
	if (node->DSW != NULL && (node->DSW != down || inverse)) {
		var = applyForce(node->DSW, body);
		if (!var) {
			forceOverNode(node->DSW, node, body, 1);
		}
	}
	if (!inverse)
		forceOverNode(node->UP, node, body, 0);
	return;
}

void init(void)
{
	int i, j;
	double vx, vy, vz, theta, dist, dist2;
	positionData = fopen("positionData.csv", "wb");
	nodes = (struct node *)malloc(sizeof(struct node) * MAX_NODES);
	roots = (struct node **)malloc(sizeof(struct node *) * MAX_NODES);
	bodies = (struct body *)malloc(sizeof(struct body) * BODIES_QUANTITY);
	for (i = 0; i < MAX_NODES; i++) {
		nodes[i].initialized = 0;
		nodes[i].UNE = NULL;
		nodes[i].UNW = NULL;
		nodes[i].USE = NULL;
		nodes[i].USW = NULL;
		nodes[i].DNE = NULL;
		nodes[i].DNW = NULL;
		nodes[i].DSE = NULL;
		nodes[i].DSW = NULL;
	}
	initializeNode(&nodes[0], NULL, -SIZE_OF_SIMULATION,
		       -SIZE_OF_SIMULATION, -SIZE_OF_SIMULATION,
		       SIZE_OF_SIMULATION, SIZE_OF_SIMULATION,
		       SIZE_OF_SIMULATION, BODIES_QUANTITY, 0);
        for(i = 0; i < BODIES_QUANTITY; i++){
		bodies[i].position(0) = i % 20 * 5;
		bodies[i].position(1) = (i / 20) % 20 * 5;
		bodies[i].position(2) = (i / 400) * 5;
		bodies[i].force(0) = 0;
		bodies[i].force(1) = 0;
		bodies[i].force(2) = 0;
		bodies[i].density = 1;
		bodies[i].pression = 1;
		bodies[i].Ss = 300;
		bodies[i].h = 20;
		theta = atan2(vy, vx) + PI / 2.0;
		bodies[i].mass = 1;
		bodies[i].speed(0) = 0;	// cos(theta) * 6.5 * 10E3 * (dist / (60E3 * LY) - 1.3) * (dist / (60E3 * LY));
		bodies[i].speed(1) = 0;	// sin(theta) * 6.5 * 10E3 * (dist / (60E3 * LY) - 1.3)* (dist / (60E3 * LY));
		bodies[i].speed(2) = 0;	//(rand() % 10000000 / 10000000.0) * 5000 - 2500;
		theta = atan2(vz, vy) + PI / 2.0;
		//bodies[i].speed(2) = sin(theta) * 0.5 * 10E3 * (dist / (60E3 * LY) - 1.1) * (dist / (60E3 * LY));
		//bodies[i].speed(1) += cos(theta) * 0.5 * 10E3 * (dist / (60E3 * LY) - 1.1) * (dist / (60E3 * LY));
		addBodyInNode(&bodies[i], &nodes[0]);
	}
}
