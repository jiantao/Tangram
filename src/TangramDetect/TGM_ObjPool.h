/*
 * =====================================================================================
 *
 *       Filename:  TGM_ObjPool.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/18/2012 11:31:24 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_OBJPOOL_H
#define  TGM_OBJPOOL_H

#include <deque>
#include <vector>

namespace Tangram
{
    template <class T> class ObjPool
    {
        public:

            ObjPool();
            ~ObjPool();

            T& Alloc(unsigned int& idx);

            inline void Return(const unsigned int& idx);

            inline unsigned int UsedSize(void) const
            {
                return pool.size() - candidates.size();
            }

            inline unsigned int Capacity(void) const
            {
                return pool.size();
            }

            inline bool IsExisted(unsigned int idx) const
            {
                return checkMap[idx];
            }

            inline T& operator[](unsigned int i)
            {
                return pool[i];
            }

            inline const T& operator[](unsigned int i) const
            {
                return pool[i];
            }

        private:

            static const unsigned int DEFAULT_SIZE = 50;

            std::deque<T> pool;

            std::deque<unsigned int> candidates;

            std::vector<bool> checkMap;
    };

    template <class T> ObjPool<T>::ObjPool()
    {
        pool.resize(DEFAULT_SIZE);
        checkMap.resize(DEFAULT_SIZE, false);
        candidates.resize(DEFAULT_SIZE);

        for (unsigned int i = 0; i != DEFAULT_SIZE; ++i)
            candidates[i] = i;
    }

    template <class T> ObjPool<T>::~ObjPool()
    {

    }

    template <class T> T& ObjPool<T>::Alloc(unsigned int& idx)
    {
        if (candidates.empty())
        {
            unsigned int oldSize = pool.size();
            unsigned int newSize = oldSize * 2;
            pool.resize(newSize);
            checkMap.resize(newSize, false);

            for (unsigned int i = oldSize - 1; i != newSize; ++i)
                candidates.push_back(i);
        }

        idx = candidates[0];
        candidates.pop_front();
        checkMap[idx] = true;
        return pool[idx];
    }

    template <class T> inline void ObjPool<T>::Return(const unsigned int& idx)
    {
        candidates.push_front(idx);
        checkMap[idx] = false;
    }
};

#endif  /*TGM_OBJPOOL_H*/
