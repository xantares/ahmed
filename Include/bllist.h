/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BLLIST_H
#define BLLIST_H

#include "blcluster.h"

extern void gen_BlSequence(blcluster*, blcluster**&);
extern void gen_upBlSequence(blcluster*, blcluster**&, unsigned&);
extern void gen_lwBlSequence(blcluster*, blcluster**&, unsigned&);
extern void gen_HilbertBlSeq(blcluster*, blcluster**&);
extern void gen_NortonBlSeq(blcluster*, blcluster**&);
extern void genBlSeqPart(blcluster*, unsigned, blcluster**&, unsigned*&,
                         unsigned (*cost_fnct)(blcluster&)=NULL);
#endif

