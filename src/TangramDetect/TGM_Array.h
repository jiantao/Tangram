/*
 * =====================================================================================
 *
 *       Filename:  TGM_Array.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/07/2012 04:23:10 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_ARRAY_H
#define  TGM_ARRAY_H

#include <iostream>
#include <cstdlib>

namespace Tangram
{
    typedef int (*CompareFunc)(const void* a, const void* b);

    template <class T> class Array
    {

        public:
            Array();
            ~Array();

            inline void Init(unsigned int capacity);
            inline void MemSet(int value);
            inline void InitToEnd(void);

            void ResizeNoCopy(unsigned int newCap);
            void Resize(unsigned int newCap);

            inline void Increment(void);
            inline void SetSize(unsigned int newSize);
            inline void Clear(void);

            inline T* GetPointer(unsigned int);
            inline T& operator[](unsigned int i);

            void Sort(CompareFunc compare);

            inline T& Last(void)
            {
                return data[size - 1];
            }

            inline T& End(void)
            {
                return data[size];
            }

        public:

            inline bool IsFull(void) const;
            inline unsigned int Size(void) const;
            inline unsigned int Capacity(void) const;

            inline int BinarySearch(const T* key, CompareFunc compare) const;

            int UpperBound(const T& key, CompareFunc compare) const;

            inline const T* GetPointer(unsigned int) const;
            inline const T& operator[](unsigned int i) const;

            inline const T& First(void) const;
            inline const T& Last(void) const;

            inline const T& End(void) const
            {
                return data[size];
            }

        private:

            Array(const Array&);
            Array& operator=(const Array&);

            T* data;
            unsigned int size;
            unsigned int capacity;
    };


    template <class T> Array<T>::Array()
    {
        data = NULL;
        size = 0;
        capacity = 0;
    }

    template <class T> Array<T>::~Array()
    {
        free(data);
        data = NULL;
    }

    template <class T> inline void Array<T>::Init(unsigned int capacity)
    {
        ResizeNoCopy(capacity);
    }

    template <class T> inline void Array<T>::InitToEnd(void)
    {
        memset(data + size, 0, sizeof(T) * (capacity - size));
    }

    template <class T> inline void Array<T>::MemSet(int value)
    {
        memset(data, value, sizeof(T) * (capacity));
    }

    template <class T> void Array<T>::ResizeNoCopy(unsigned int newCap)
    {
        if (capacity < newCap)
        {
            free(data);
            data = (T*) calloc(sizeof(T), newCap);
            if (data == NULL)
            {
                std::cerr << "ERROR: Not enough memory for the new element in the array.\n";
                exit(1);
            }

            capacity = newCap;
        }

        size = 0;
    }

    template <class T> void Array<T>::Resize(unsigned int newCap)
    {
        data = (T*) realloc(data, sizeof(T) * newCap);
        if (data == NULL)
        {
            std::cerr << "ERROR: Not enough memory for the new element in the array.\n";
            exit(1);
        }

        capacity = newCap;

        if (size > newCap)
            size = newCap;
    }

    template <class T> inline void Array<T>::Increment(void)
    {
        ++size;
    }

    template <class T> inline void Array<T>::SetSize(unsigned int newSize)
    {
        size = newSize;
    }

    template <class T> inline void Array<T>::Clear(void)
    {
        size = 0;
    }

    template <class T> T* Array<T>::GetPointer(unsigned int i)
    {
        return data + i;
    }

    template <class T> T& Array<T>::operator[](unsigned int i)
    {
        return data[i];
    }


    template <class T> void Array<T>::Sort(CompareFunc compare)
    {
        qsort(data, size, sizeof(T), compare);
    }

    template <class T> inline bool Array<T>::IsFull(void) const 
    {
        return (size == capacity);
    }

    template <class T> inline unsigned int Array<T>::Size(void) const
    {
        return size;
    }

    template <class T> inline unsigned int Array<T>::Capacity(void) const
    {
        return capacity;
    }

    template <class T> inline int Array<T>::BinarySearch(const T* pKey, CompareFunc compare) const
    {
        T* pValue = ((T*) bsearch(pKey, data, size, sizeof(T), compare));
        if (pValue != NULL)
            return (pValue - data);
        else
            return -1;
    }

    template <class T> int Array<T>::UpperBound(const T& key, CompareFunc compare) const
    {
        int low = 0;
        int high = size - 1;

        while (low < high)
        {
            int mid = (low + high) / 2;
            if (compare(&key, GetPointer(mid)) <= 0)
                high = mid;
            else
                low = mid + 1;
        }

        if (compare(&key, GetPointer(low)) <= 0)
            return low;
        else
            return -1;
    }

    template <class T> inline const T* Array<T>::GetPointer(unsigned int i) const
    {
        return data + i;
    }

    template <class T> inline const T& Array<T>::operator[](unsigned int i) const
    {
        return data[i];
    }

    template <class T> inline const T& Array<T>::First(void) const
    {
        return data[0];
    }

    template <class T> inline const T& Array<T>::Last(void) const
    {
        return data[size - 1];
    }
};


#endif  /*TGM_ARRAY_H*/
