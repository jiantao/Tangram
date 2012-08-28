/*
 * =====================================================================================
 *
 *       Filename:  TGM_Cluster.h
 *
 *    Description:  2D Nearest Neighbours Clustering
 *
 *        Version:  1.0
 *        Created:  05/11/2012 10:31:53 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 * =====================================================================================
 */

#ifndef  TGM_CLUSTER_H
#define  TGM_CLUSTER_H

#include "TGM_Array.h"
#include "TGM_PairAttrbt.h"

namespace Tangram
{
    struct ClusterElmnt
    {
        double mean[2];              // mean of the first and the second attribute

        double std[2];               // standard deviatioin of the fisrt and the second attribute

        double min[2];               // minimum values of the first and the second attribute

        double max[2];               // max values of the first attribute and the second attribute

        unsigned int numReadPair;    // number of read pairs

        unsigned int startIndex;     // start index of the cluster elements
    };

    class Cluster
    {
        public:
            Cluster();
            ~Cluster();

            // initialize the clustering elements
            void Init(const Array<PairAttrbt>* pPairAttabts, unsigned int minClusterSize, const double* minStd);

            // make the clusters
            void Make(void);

            // print the cluster elements
            void Print(int32_t refID, const char* specialID) const;

            inline const Array<ClusterElmnt>* GetCluterElmnts(void) const
            {
                return &elements;
            }

            inline unsigned int Size(void) const
            {
                return elements.Size();
            }

            inline unsigned int GetActualNumElmnts(void) const
            {
                return actualNumElmnts;
            }

            inline const Array<unsigned int>& GetNextArray(void) const
            {
                return next;
            }

            // connect a member to its center
            void Connect(unsigned int centerIndex, unsigned int memberIndex);

        private:

            // build the basic elements for clustering
            void Build(void);

            // create the final clusters
            void Finalize(void);

            // add a new cluster element
            unsigned int AddElmnt(unsigned int attrbtIdx);

            // update the existing cluster element
            void UpdateElmnt(unsigned int elmntIdx, unsigned int attrbtIdx);

            // clean the clusters
            void Clean(void);

        public:

            // a pointer to the pair attribute arrays
            const Array<PairAttrbt>* pPairAttrbts;

            // helper array to build a circled linked list
            Array<unsigned int> next;

        private:

            // array of cluster elements
            Array<ClusterElmnt> elements;

            // a map from one point to its center point
            Array<int> map;

            // number of neighours around each point
            Array<unsigned int> count;

            // actual number of cluster elements
            unsigned int actualNumElmnts;

            // minimum number of members in a cluster
            unsigned int minClusterSize;

            // minimum standard deviation of a cluster
            double minStd[2];
    };
};

#endif  /*TGM_CLUSTER_H*/
