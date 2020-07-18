#include "memory.h"

using namespace std;


float ** allocateTwoDim_float(unsigned row, unsigned col)
{
    float ** ptr = new float*[row];
    for(int i = 0; i < row; i++)
    {
        ptr[i] = new float[col];
        for (unsigned j(0) ; j < col; ++j)
        {
            ptr[i][j] = float(0);
        }

    }
    return ptr;
}

/****************************************************************************************************************************/
/****************************************************************************************************************************/



short int ** allocateTwoDim_short(unsigned int row, unsigned int col)
{
    short int ** ptr = new short int*[row];
    for(int i = 0; i < row; i++)
    {
        ptr[i] = new short int[col];
        for (unsigned j(0) ; j < col; ++j)
        {
            ptr[i][j] = float(0);
        }
    }
    return ptr;
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

void destroyTwoDim_short(short ** ptr, unsigned row, unsigned col)
{
    for(int i = 0; i < row; i++)
    {
        delete [] ptr[i];
    }
    delete [] ptr;
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

void destroyTwoDim_float(float ** ptr, unsigned row, unsigned col)
{
    for(int i = 0; i < row; i++)
    {
        delete [] ptr[i];
    }
    delete [] ptr;
}


